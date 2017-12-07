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

#include <klib/extern.h>
#include <klib/hashtable.h>
#include <klib/rc.h>
#include <klib/sort.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef memcpy

#ifdef __cplusplus
extern "C" {
#endif

static const uint64_t BUCKET_VALID = (uint64_t)1 << 63;

LIB_EXPORT rc_t CC KHashTableInit(KHashTable* self, size_t key_size,
                                  size_t value_size, size_t initial_buckets,
                                  double max_load_factor, bool key_cstr)
{
    if (self == NULL)
        return RC(rcCont, rcHashtable, rcConstructing, rcParam, rcInvalid);

    if (max_load_factor < 0)
        return RC(rcCont, rcHashtable, rcConstructing, rcParam, rcInvalid);

    if (max_load_factor == 0.0)
        self->max_load_factor = 0.6;
    else
        self->max_load_factor = max_load_factor;

    if (key_size == 0)
        return RC(rcCont, rcHashtable, rcConstructing, rcParam, rcInvalid);

    if (key_cstr && key_size != 8)
        return RC(rcCont, rcHashtable, rcConstructing, rcParam, rcInvalid);

    if (initial_buckets == 0) {
        initial_buckets = 16;
    } else {
        uint32_t lg2 = (64 - __builtin_clzll(initial_buckets | 1));
        initial_buckets = 1 << lg2;
    }

    self->key_size = key_size;
    self->value_size = value_size;
    self->num_buckets = initial_buckets;
    self->mask = initial_buckets - 1;
    self->count = 0;
    self->key_cstr = key_cstr;

    self->bucket_size = 8 + key_size + value_size;
    // keyhash (u64) + key + value, MSB of byckethash initially set to 0 for
    // invalid
    self->buckets = calloc(1, self->num_buckets * self->bucket_size);
    if (!self->buckets)
        return RC(rcCont, rcHashtable, rcConstructing, rcMemory, rcExhausted);

    return 0;
}

LIB_EXPORT void CC
KHashTableWhack(KHashTable* self, void(CC* keywhack)(void* item, void* data),
                void(CC* valuewhack)(void* item, void* data), void* data)
{
    if (self == NULL) return;

    // TODO: for ...
    if (keywhack != NULL) {
    }
    if (valuewhack != NULL) {
    }
    self->key_size = 0;
    self->value_size = 0;
    self->num_buckets = 0;
    self->count = 0;
    self->mask = 0;
    free(self->buckets);
    self->buckets = NULL;
    memset(self, 0, sizeof(KHashTable));
}

LIB_EXPORT size_t CC KHashTableCount(const KHashTable* self)
{
    if (self != NULL)
        return self->count;
    else
        return 0;
}

static rc_t resize(KHashTable* self)
{
    void* old_buckets = self->buckets;
    size_t old_num_buckets = self->num_buckets;
    self->num_buckets *= 2;
    void* new_buckets = calloc(1, self->num_buckets * self->bucket_size);
    if (!new_buckets)
        return RC(rcCont, rcHashtable, rcInserting, rcMemory, rcExhausted);

    self->buckets = new_buckets;
    self->mask = (self->mask << 1) | 1;
    self->count = 0;

    for (size_t bucket = 0; bucket != old_num_buckets; ++bucket) {
        const char* bucketptr
            = (char*)old_buckets + (bucket * self->bucket_size);
        const char* hashptr = bucketptr;
        const char* keyptr = bucketptr + 8;
        const char* valueptr = bucketptr + 8 + self->key_size;
        uint64_t buckethash;
        memcpy(&buckethash, hashptr, 8);
        if (buckethash & BUCKET_VALID) {
            KHashTableAdd(self, keyptr, buckethash, valueptr);
        }
    }

    free(old_buckets);
    return 0;
}

/*
Since no C++ templates, buckets is going to be variable sized struct:
{
    u64 keyhash; // we use MSB=1 to indicate validity
    u8[key_size] key;
    u8[value_size] value; // Will not be present for sets
}
*/

LIB_EXPORT bool CC KHashTableFind(const KHashTable* self, const void* key,
                                  uint64_t keyhash, void* value)

{
    if (self == NULL || self->buckets == NULL) return false;

    keyhash |= BUCKET_VALID;
    uint64_t bucket = keyhash;
    while (1) {
        bucket &= self->mask;
        const char* bucketptr
            = (char*)self->buckets + (bucket * self->bucket_size);
        const char* hashptr = bucketptr;
        const char* keyptr = bucketptr + 8;
        const char* valueptr = bucketptr + 8 + self->key_size;
        uint64_t buckethash;
        memcpy(&buckethash, hashptr, 8);
        if (!(buckethash & BUCKET_VALID)) // reached invalid bucket
        {
            return false;
        }

        if (buckethash == keyhash)
        // hash hit, odds are very low (2^-63) that this is an actual miss,
        // but we have to check.
        {
            bool found;

            if (self->key_cstr) {
                char* p;
                memcpy(&p, keyptr, 8);
                found = (strcmp(p, key) == 0);
            } else
                found = (memcmp(keyptr, key, self->key_size) == 0);

            if (found) {
                if (value && self->value_size)
                    memcpy(value, valueptr, self->value_size);

                return true;
            }
        }
        ++bucket;
    }
}

LIB_EXPORT rc_t CC KHashTableAdd(KHashTable* self, const void* key,
                                 uint64_t keyhash, const void* value)
{
    if (self == NULL || self->buckets == NULL)
        return RC(rcCont, rcHashtable, rcInserting, rcParam, rcInvalid);

    keyhash |= BUCKET_VALID;
    uint64_t bucket = keyhash;
    while (1) {
        bucket &= self->mask;
        char* bucketptr = (char*)self->buckets + (bucket * self->bucket_size);
        char* hashptr = bucketptr;
        char* keyptr = bucketptr + 8;
        char* valueptr = bucketptr + 8 + self->key_size;
        uint64_t buckethash;
        memcpy(&buckethash, bucketptr, 8);
        if (!(buckethash & BUCKET_VALID)) // reached invalid bucket
        {
            memcpy(hashptr, &keyhash, 8);

            if (self->key_cstr)
                memcpy(keyptr, &key, 8);
            else
                memcpy(keyptr, key, self->key_size);

            if (self->value_size) {
                memcpy(valueptr, value, self->value_size);
            }

            ++self->count;

            if (KHashTableGetLoadFactor(self) > self->max_load_factor)
                return resize(self);

            return 0;
        }

        if (buckethash == keyhash) // hash hit
        {
            bool found;
            if (self->key_cstr) {
                char* p;
                memcpy(&p, keyptr, 8);
                found = (strcmp(p, key) == 0);
            } else {
                found = (memcmp(keyptr, key, self->key_size) == 0);
            }

            if (found) {
                // replacement
                if (self->value_size)
                    memcpy(valueptr, value, self->value_size);
                return 0;
            }
        }

        ++bucket;
    }
}

// Return current load factor (# buckets / # items)
LIB_EXPORT double CC KHashTableGetLoadFactor(const KHashTable* self)
{
    if (self == NULL) return 0.0;

    double load_factor = (double)self->count / (double)self->num_buckets;
    return load_factor;
}

// TODO: Iterator, can become invalid after any insert

// Fast Hash function
// Inner core from Google's FarmHash.
// Removed
//  * 32-bit support
//  * SSE4.1/AVX requirement
//  * big-endian
//  * long string optimizations
//  * STL
//  * Seeding

// Some primes between 2^63 and 2^64 for various uses.
static const uint64_t k0 = 0xc3a5c85c97cb3127ULL;
static const uint64_t k1 = 0xb492b66fbe98f273ULL;
static const uint64_t k2 = 0x9ae16a3b2f90404fULL;

#if !defined(uint128_t)
#define uint128_t __uint128_t
#endif

uint64_t ShiftMix(uint64_t val) { return val ^ (val >> 47); }

static inline uint64_t Uint128Low64(const uint128_t x)
{
    return (uint64_t)(x);
}
static inline uint64_t Uint128High64(const uint128_t x)
{
    return (uint64_t)(x >> 64);
}
static inline uint128_t Uint128(uint64_t lo, uint64_t hi)
{
    return lo + (((uint128_t)hi) << 64);
}

/*
inline uint64_t Hash128to64(uint128_t x)
{
    // Murmur-inspired hashing.
    const uint64_t kMul = 0x9ddfea08eb382d69ULL;
    uint64_t a = (Uint128Low64(x) ^ Uint128High64(x)) * kMul;
    a ^= (a >> 47);
    uint64_t b = (Uint128High64(x) ^ a) * kMul;
    b ^= (b >> 47);
    b *= kMul;
    return b;
}
*/

static inline uint64_t HashLen16(uint64_t u, uint64_t v, uint64_t mul)
{
    // Murmur-inspired hashing.
    uint64_t a = (u ^ v) * mul;
    a ^= (a >> 47);
    uint64_t b = (v ^ a) * mul;
    b ^= (b >> 47);
    b *= mul;
    return b;
}

static inline uint64_t Fetch64(const char* p)
{
    uint64_t result;
    memcpy(&result, p, sizeof(result));
    return result;
}

static inline uint32_t Fetch32(const char* p)
{
    uint32_t result;
    memcpy(&result, p, sizeof(result));
    return result;
}

static inline uint64_t Rotate64(uint64_t val, int shift)
{
    // Avoid shifting by 64: doing so yields an undefined result.
    return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
}

static inline uint64_t HashLen0to16(const char* s, size_t len)
{
    if (len >= 8) {
        uint64_t mul = k2 + len * 2;
        uint64_t a = Fetch64(s) + k2;
        uint64_t b = Fetch64(s + len - 8);
        uint64_t c = Rotate64(b, 37) * mul + a;
        uint64_t d = (Rotate64(a, 25) + b) * mul;
        return HashLen16(c, d, mul);
    }
    if (len >= 4) {
        uint64_t mul = k2 + len * 2;
        uint64_t a = Fetch32(s);
        return HashLen16(len + (a << 3), Fetch32(s + len - 4), mul);
    }
    if (len > 0) {
        uint8_t a = s[0];
        uint8_t b = s[len >> 1];
        uint8_t c = s[len - 1];
        uint32_t y = (uint32_t)(a) + ((uint32_t)(b) << 8);
        uint32_t z = len + ((uint32_t)(c) << 2);
        return ShiftMix(y * k2 ^ z * k0) * k2;
    }
    return k2;
}

static inline uint64_t HashLen17to32(const char* s, size_t len)
{
    uint64_t mul = k2 + len * 2;
    uint64_t a = Fetch64(s) * k1;
    uint64_t b = Fetch64(s + 8);
    uint64_t c = Fetch64(s + len - 8) * mul;
    uint64_t d = Fetch64(s + len - 16) * k2;
    return HashLen16(Rotate64(a + b, 43) + Rotate64(c, 30) + d,
                     a + Rotate64(b + k2, 18) + c, mul);
}

static inline uint64_t HashLen33to64(const char* s, size_t len)
{
    uint64_t mul = k2 + len * 2;
    uint64_t a = Fetch64(s) * k2;
    uint64_t b = Fetch64(s + 8);
    uint64_t c = Fetch64(s + len - 8) * mul;
    uint64_t d = Fetch64(s + len - 16) * k2;
    uint64_t y = Rotate64(a + b, 43) + Rotate64(c, 30) + d;
    uint64_t z = HashLen16(y, a + Rotate64(b + k2, 18) + c, mul);
    uint64_t e = Fetch64(s + 16) * mul;
    uint64_t f = Fetch64(s + 24);
    uint64_t g = (y + Fetch64(s + len - 32)) * mul;
    uint64_t h = (z + Fetch64(s + len - 24)) * mul;
    return HashLen16(Rotate64(e + f, 43) + Rotate64(g, 30) + h,
                     e + Rotate64(f + a, 18) + g, mul);
}

LIB_EXPORT uint64_t CC KHash(const char* s, size_t len)
{
    if (len <= 32) {
        if (len <= 16) {
            return HashLen0to16(s, len);
        } else {
            return HashLen17to32(s, len);
        }
    } else if (len <= 64) {
        return HashLen33to64(s, len);
    } else {
        uint64_t h = 0;
        while (len >= 64) {
            h += HashLen33to64(s, 64);
            s += 64;
            len -= 64;
        }
        return h + KHash(s, len);
    }
}

#ifdef __cplusplus
}
#endif
