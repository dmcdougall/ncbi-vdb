/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/
#include <arch-impl.h>
#include <assert.h>
#include <atomic.h>
#include <kfs/mmap.h>
#include <klib/hashfile.h>
#include <klib/rc.h>
#include <klib/vector.h>
#include <kproc/lock.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef memcpy

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char u8;
static u8* BUCKET_INVALID = (u8*)0; /* Must be 0, calloc fills */
static u8* BUCKET_INVISIBLE = (u8*)1;
#define NUM_SEGMENTS 64

#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#if defined(__GNUC__)
#define STORE_FENCE __sync_synchronize()
#elif defined(_MSC_VER)
#define STORE_FENCE _ReadWriteBarrier()
#else
/* C11's stdatomic */
#define STORE_FENCE atomic_thread_fence(memory_order_release)
#endif

typedef struct HKV
{
    uint64_t hash;
    size_t key_size;
    size_t value_size;
    void* key;
    void* value;
} HKV;

typedef struct Hashtable
{
    u8** table;
    size_t table_sz;
} Hashtable;

typedef struct Segment
{
    Hashtable* hashtable; /* Will be switched atomically so readers can be
                             lock-free */
    size_t load;          /* Including invisible entries */
    KLock* seglock;       /* TODO, use adaptive spinlocks for linux */
    u8* alloc_base;
    size_t alloc_remain;
    bool in_rehash;
} Segment;

struct KHashFile
{
    KFile* file;
    atomic64_t count;
    int64_t iterator;
    Segment segments[NUM_SEGMENTS];
    KLock* alloc_lock; /* protects below */
    Vector allocs;
};

/* Hash/Key/Values are encoded as:
 * 1 byte key length, if key length >255 store 255 and actual size below
 * 1 byte value length, if value length >255 store 255 and actual size below
 * 8 bytes hash - 6 enough for a 256 trillion entry set requiring 5PB.
 * [ 8 bytes ] if key length >255, store length here
 * [ 8 bytes ] if value length >255, store length here
 * key
 * value
 */
static size_t hkv_space(const HKV* in)
{
    size_t out = 0;
    out += 8;                           /* hash */
    out += 2;                           /* key_size/value_size */
    if (in->key_size > 254) out += 8;   /* optional key_size */
    if (in->value_size > 254) out += 8; /* optional value_size */
    out += in->key_size;                /* key */
    out += in->value_size;              /* value */
    return out;
}

static void hkv_decode(const u8* in, HKV* out)
{
    u8* p = (u8*)in;
    memcpy(&out->hash, p, 8);
    p += 8;
    out->key_size = *p;
    p += 1;
    out->value_size = *p;
    p += 1;
    if (out->key_size == 255) {
        memcpy(&out->key_size, p, 8);
        p += 8;
    }
    if (out->value_size == 255) {
        memcpy(&out->value_size, p, 8);
        p += 8;
    }
    out->key = p;
    p += out->key_size;
    out->value = p;
    p += out->value_size;
}

static void hkv_encode(const HKV* in, u8* out)
{
    u8* p = (u8*)out;
    memcpy(p, &in->hash, 8);
    p += 8;
    if (in->key_size <= 254)
        *p = (u8)in->key_size;
    else
        *p = 255;
    ++p;
    if (in->value_size <= 254)
        *p = (u8)in->value_size;
    else
        *p = 255;
    ++p;
    if (in->key_size > 255) {
        memcpy(out, &in->key_size, 8);
        p += 8;
    }
    if (in->value_size > 255) {
        memcpy(p, &in->value_size, 8);
        p += 8;
    }
    memcpy(p, in->key, in->key_size);
    p += in->key_size;
    memcpy(p, in->value, in->value_size);
    p += in->value_size;
}

#if 0
static void dump(const KHashFile * self)
{
    assert(self != NULL);
    if (self->num_buckets > 40) return;
    fprintf(stderr, "-- table has %ld/%ld\n", self->count, self->load);
    for (size_t bucket = 0; bucket != self->num_buckets; bucket++) {
        const char* bucketptr = (char*)self->buckets;
        const char* hashptr = bucketptr;
        const char* keyptr = bucketptr + 8;
        uint64_t buckethash=0;
        fprintf(stderr, "   bucket %03ld hash %lx", bucket, buckethash);
        if (buckethash & BUCKET_VALID) fprintf(stderr, " val");
        if (buckethash & BUCKET_VISIBLE) fprintf(stderr, " vis");
        /*fprintf(stderr, " key=%ld", key); */
        fprintf(stderr, "\n");
    }
}
#endif

/* Single, locked allocator shared between all segments */
static void* map_calloc(KHashFile* self, size_t size, size_t nmemb)
{
    rc_t rc;
    void* block = NULL;

    if (self == NULL || size == 0) return NULL;

    rc = KLockAcquire(self->alloc_lock);
    if (rc) return NULL;

    /* Round up to 4K page size */
    size *= nmemb;
    size = ((size + 4095) / 4096) * 4096;

    if (self->file != NULL) {
        KMMap* mm;

        uint64_t filesize;
        rc = KFileSize(self->file, &filesize);
        if (rc) {
            KLockUnlock(self->alloc_lock);
            return NULL;
        }
        rc = KFileSetSize(self->file, filesize + size);
        if (rc) {
            KLockUnlock(self->alloc_lock);
            return NULL;
        }
        rc = KMMapMakeRgnUpdate(&mm, self->file, filesize, size);
        if (rc) {
            KLockUnlock(self->alloc_lock);
            return NULL;
        }
        rc = KMMapAddrUpdate(mm, &block);
        if (rc) {
            KLockUnlock(self->alloc_lock);
            return NULL;
        }
    } else {
        block = calloc(1, size);
        VectorAppend(&self->allocs, NULL, block);
    }

    rc = KLockUnlock(self->alloc_lock);
    if (rc) return NULL;

    return block;
}

/* Per segment allocator */
static void* seg_alloc(KHashFile* self, size_t segment, size_t size)
{
    assert(self != NULL);
    assert(segment < NUM_SEGMENTS);
    assert(size > 0);

    Segment* seg = &self->segments[segment];
    if (size > seg->alloc_remain) {
        size_t req = MAX(size, 65536);
        seg->alloc_base = map_calloc(self, req, 1);
        seg->alloc_remain = req;
    }

    void* r = seg->alloc_base;
    seg->alloc_base += size;
    seg->alloc_remain -= size;

    return r;
}

static rc_t rehash(KHashFile* self, size_t segment, size_t capacity)
{
    assert(self != NULL);
    assert(segment < NUM_SEGMENTS);

    /* Assume segment is locked by caller */

    self->iterator = -1; /* Invalidate any current iterators */

    Segment* seg = &self->segments[segment];

    Hashtable* old_hashtable = seg->hashtable;

    /* Don't allow shrinking */
    if (old_hashtable) capacity = MAX(capacity, old_hashtable->table_sz);

    uint64_t lg2 = (uint64_t)uint64_msbit(capacity | 1);
    capacity = 1ULL << lg2;

    Hashtable* hashtable = (Hashtable*)map_calloc(self, 1, sizeof(Hashtable));
    if (!hashtable)
        return RC(rcCont, rcHashtable, rcInserting, rcMemory, rcExhausted);

    hashtable->table_sz = capacity;
    hashtable->table = map_calloc(self, capacity + 1, sizeof(u8*));
    if (!hashtable->table) {
        return RC(rcCont, rcHashtable, rcInserting, rcMemory, rcExhausted);
    }

    if (old_hashtable) {
        seg->load = 0;
        seg->in_rehash = true;
        u8** table = old_hashtable->table;
        for (size_t bucket = 0; bucket != old_hashtable->table_sz; ++bucket) {
            u8* kv = table[bucket];
            if (kv != BUCKET_INVALID && kv != BUCKET_INVISIBLE) {
                HKV bkv;
                hkv_decode(kv, &bkv);
                KHashFileAdd(self, bkv.key, bkv.key_size, bkv.hash, bkv.value,
                             bkv.value_size);
            }
        }
        seg->in_rehash = false;
    }
    STORE_FENCE;
    seg->hashtable = hashtable;
    free(old_hashtable);

    return 0;
}

LIB_EXPORT rc_t KHashFileMake(KHashFile** self, KFile* hashfile)
{
    if (self == NULL)
        return RC(rcCont, rcHashtable, rcConstructing, rcParam, rcInvalid);

    rc_t rc;
    *self = NULL;

    KHashFile* kht = (KHashFile*)malloc(sizeof(KHashFile));
    if (kht == NULL)
        return RC(rcCont, rcHashtable, rcConstructing, rcMemory, rcExhausted);

    kht->file = hashfile;
    atomic64_set(&kht->count, 0);
    VectorInit(&kht->allocs, 0, 0);
    rc = KLockMake(&kht->alloc_lock);
    if (rc) return rc;
    for (size_t i = 0; i != NUM_SEGMENTS; ++i) {
        kht->segments[i].hashtable = NULL;
        kht->segments[i].load = 0;
        rc = KLockMake(&kht->segments[i].seglock);
        if (rc) return rc;
        kht->segments[i].alloc_base = NULL;
        kht->segments[i].alloc_remain = 0;
        rc = rehash(kht, i, 4);

        if (rc) {
            free(kht);
            return rc;
        }
    }

    *self = kht;

    return 0;
}

LIB_EXPORT void KHashFileDispose(KHashFile* self)
{
    if (self == NULL) return;

    /* TODO: Stop the world? */
    /*    KLockWhack(&self->alloc_lock);  */
    for (size_t i = 0; i != NUM_SEGMENTS; ++i) {
        self->segments[i].alloc_base = NULL;
        self->segments[i].alloc_remain = 0;
        /*        KLockWhack(&self->segments[i].seglock); */
    }

    for (uint32_t i = 0; i != VectorLength(&self->allocs); ++i) {
        void* alloc = VectorGet(&self->allocs, i);
        free(alloc);
    }
    /* TODO: munmap? */
    VectorWhack(&self->allocs, NULL, NULL);

    memset(self, 0, sizeof(KHashFile));
    self->iterator = -1;
    free(self);
}

LIB_EXPORT size_t KHashFileCount(const KHashFile* self)
{
    if (self != NULL)
        return atomic64_read(&self->count);
    else
        return 0;
}

LIB_EXPORT bool KHashFileFind(const KHashFile* self, const void* key,
                              const size_t key_size, const uint64_t keyhash,
                              void* value, size_t* value_size)
{
    if (self == NULL)
        return RC(rcCont, rcHashtable, rcInserting, rcParam, rcInvalid);

    size_t triangle = 0;
    uint64_t bucket = keyhash;
    int segment = bucket & (NUM_SEGMENTS - 1);
    const Segment* seg = &self->segments[segment];
    Hashtable* hashtable = seg->hashtable;
    /* TODO: Check asm that hashtable isn't reloaded */
    STORE_FENCE;
    u8** table = hashtable->table;
    const uint64_t mask = hashtable->table_sz - 1;

    while (1) {
        bucket &= mask;
        u8* kv = table[bucket];
        if (kv == BUCKET_INVALID) return false;

        if (kv != BUCKET_INVISIBLE) {
            uint64_t buckethash = *(uint64_t*)kv;

            if (buckethash == keyhash) /* hash hit */
            {
                HKV bkv;
                hkv_decode(kv, &bkv);
                if (bkv.key_size == key_size
                    && (memcmp(key, bkv.key, key_size) == 0)) {
                    if (value) memcpy(value, bkv.value, bkv.value_size);
                    if (value_size) *value_size = bkv.value_size;
                    return true;
                }
            }
        }

        /* To improve lookups when hash function has poor distribution, we use
         * quadratic probing with a triangle sequence: 0,1,3,6,10...
         * This will allow complete coverage on a % 2^N hash table.
         */
        ++triangle;
        bucket += (triangle * (triangle + 1) / 2);
    }
}

LIB_EXPORT rc_t KHashFileAdd(KHashFile* self, const void* key,
                             const size_t key_size, const uint64_t keyhash,
                             const void* value, const size_t value_size)
{
    if (self == NULL)
        return RC(rcCont, rcHashtable, rcInserting, rcParam, rcInvalid);

    size_t triangle = 0;
    uint64_t bucket = keyhash;
    int segment = bucket & (NUM_SEGMENTS - 1);
    Segment* seg = &self->segments[segment];
    rc_t rc;

    if (!seg->in_rehash) {
        /* Lock segment, TODO:adaptive spin lock best */
        rc = KLockAcquire(seg->seglock);
        if (rc) return rc;
    }

    Hashtable* hashtable = seg->hashtable;
    STORE_FENCE;
    u8** table = hashtable->table;
    const uint64_t mask = hashtable->table_sz - 1;
    HKV hkv;
    hkv.hash = keyhash;
    hkv.key_size = key_size;
    hkv.value_size = value_size;
    hkv.key = (void*)key;
    hkv.value = (void*)value;
    size_t kvsize = hkv_space(&hkv);

    while (1) {
        bucket &= mask;
        u8* kv = table[bucket];
        if (kv == BUCKET_INVALID) {
            void* buf = seg_alloc(self, segment, kvsize);
            hkv_encode(&hkv, buf);

            /* Intel 64 and IA-32 Architectures Software Developers Manual
             * Volume 3 (System Programming Guide)
             * 8.2.2 Memory Ordering in P6 and More Recent Processor Families:
             *  * Writes to memory are not reordered with other writes, ...
             *
             * So as long as we update buckets last, lock-free readers
             * (aka the HashFileFinds will never see a partial key/value.
             *
             * We do need to ensure the _compiler_ doesn't reorder the stores
             * however, which is what the fence is for.
             */
            STORE_FENCE;
            table[bucket] = buf;
            seg->load++;
            if (!seg->in_rehash) {
                atomic64_inc(&self->count);
                KLockUnlock(seg->seglock);
            }

            return 0;
        }

        if (kv != BUCKET_INVISIBLE) {
            uint64_t buckethash = *(uint64_t*)kv;

            if (buckethash == keyhash) /* hash hit */
            {
                HKV bkv;
                hkv_decode(kv, &bkv);
                if (bkv.key_size == key_size
                    && (memcmp(key, bkv.key, key_size) == 0)) {
                    /* replacement */
                    if (value_size) {
                        void* buf = seg_alloc(self, segment, kvsize);
                        hkv_encode(&hkv, buf);

                        STORE_FENCE;
                        table[bucket] = buf;
                        KLockUnlock(seg->seglock);
                    }
                    return 0;
                }
            }
        }

        ++triangle;
        bucket += (triangle * (triangle + 1) / 2);
    }
}

LIB_EXPORT bool KHashFileDelete(KHashFile* self, const void* key,
                                const size_t key_size, uint64_t keyhash)
{
    if (self == NULL)
        return RC(rcCont, rcHashtable, rcInserting, rcParam, rcInvalid);

    size_t triangle = 0;
    /* Lock segment, TODO:adaptive spin lock best */
    uint64_t bucket = keyhash;
    int segment = bucket & (NUM_SEGMENTS - 1);
    Segment* seg = &self->segments[segment];

    rc_t rc = KLockAcquire(seg->seglock);
    if (rc) return false;

    Hashtable* hashtable = seg->hashtable;
    STORE_FENCE;
    u8** table = hashtable->table;
    const uint64_t mask = hashtable->table_sz - 1;

    while (1) {
        bucket &= mask;
        u8* kv = table[bucket];
        if (kv == BUCKET_INVALID) {
            KLockUnlock(seg->seglock);

            return false;
        }

        if (kv != BUCKET_INVISIBLE) {
            uint64_t buckethash = *(uint64_t*)kv;

            if (buckethash == keyhash) /* hash hit */
            {
                HKV bkv;
                hkv_decode(kv, &bkv);
                if (bkv.key_size == key_size
                    && (memcmp(key, bkv.key, key_size) == 0)) {
                    table[bucket] = BUCKET_INVISIBLE;
                    KLockUnlock(seg->seglock);
                    return true;
                }
            }
        }

        ++triangle;
        bucket += (triangle * (triangle + 1) / 2);
    }
}

LIB_EXPORT rc_t KHashFileReserve(KHashFile* self, size_t capacity)
{
    /*return rehash(self, capacity); */
    return 0;
}

LIB_EXPORT void KHashFileIteratorMake(KHashFile* self)
{
    if (self != NULL) self->iterator = 0;
}

LIB_EXPORT bool KHashFileIteratorNext(KHashFile* self, void* key,
                                      size_t* key_size, void* value,
                                      size_t* value_size)
{
    if (self == NULL || self->iterator == -1) return false;
    /*
        char* buckets = (char*)self->buckets;

        const size_t key_size = self->key_size;
        const size_t value_size = self->value_size;

        while (1) {
            if ((size_t)self->iterator >= self->num_buckets) {
                self->iterator = -1;
                return false;
            }

            char* bucketptr = buckets + (self->iterator * sizeof(uint64_t *));
            char* hashptr = bucketptr;
            const char* keyptr = bucketptr + 8;
            const char* valueptr = bucketptr + 8 + self->key_size;
            uint64_t buckethash;
            memcpy(&buckethash, hashptr, 8);

            ++self->iterator;

            if ((buckethash & BUCKET_VALID) && (buckethash & BUCKET_VISIBLE))
       {
                memcpy(key, keyptr, key_size);
                if (value && value_size) memcpy(value, valueptr, value_size);
                return true;
            }
        }
        */
    return false;
}

#ifdef __cplusplus
}
#endif
