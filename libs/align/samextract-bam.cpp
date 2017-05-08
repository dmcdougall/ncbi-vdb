/* ===========================================================================
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

#include <kfs/file.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <klib/text.h>
#include <kproc/queue.h>
#include <kproc/thread.hpp>
#include <kproc/timeout.h>
#include <kproc/lock.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strtol.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include "samextract.h"
#include "samextract-tokens.h"
#include <align/samextract-lib.h>

typedef enum BGZF_state
{
    empty,
    compressed,
    uncompressed
} BGZF_state;
typedef struct BGZF_s
{
    KLock* lock;
    Bytef in[READBUF_SZ + 1024];
    Bytef out[READBUF_SZ + 1024];
    uInt insize;
    uInt outsize;
    BGZF_state state;
} BGZF;

class BGZFview
{
  public:
    BGZFview() : bgzf(NULL), cur(NULL){};

    ~BGZFview()
    {
        releasebuf();
        bgzf = NULL;
        cur = NULL;
    }

// TODO: Handle _MSC_VER
// TODO: Move to a USE_ include that defines features:
// __HAS_RVALUE_REFERENCES, ...
// __HAS_METHOD_DELETE, ...
// __has_cpp_attribute
// __CPPVER=98,11,14,..
#ifndef __cplusplus
    // C++98
    BGZFview(const BGZFview&);
#elif __cplusplus <= 199711L
    // C++98
    BGZFview(const BGZFview&);
#elif __cplusplus >= 201103L
    // C++11, C++14
    BGZFview(const BGZFview&) = delete;             // No copy ctor
    BGZFview& operator=(const BGZFview&) = delete;  // No assignment
    BGZFview(const BGZFview&&) = delete;            // No move ctor
    BGZFview& operator=(const BGZFview&&) = delete; // No move assignment
#endif

  private:
    void releasebuf()
    {
        if (bgzf) {
            DBG("releasing");
            KLockUnlock(bgzf->lock);
            KLockRelease(bgzf->lock);
            bgzf->state = empty;
            free(bgzf);
            bgzf = NULL;
        }
    }

    bool getnextBGZF(KQueue* que)
    {
        releasebuf();

        struct timeout_t tm;
        TimeoutInit(&tm, 1000); // 1 second

        while (true)
        {
            void* where = NULL;
            DBG("popping");
            rc_t rc = KQueuePop(que, &where, &tm);
            DBG("popped");
            if (rc == 0) {
                bgzf = (BGZF*)where;
                DBG("Acquiring parser lock");
                KLockAcquire(bgzf->lock); // ready for parsing?
                DBG("Acquired parser lock");

                if (bgzf->state != uncompressed) {
                    ERR("\t\tParser queue bad state");
                    RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                       rcUnexpected);
                    return false;
                }

                DBG("\t\tParser BGZF %p size %u", bgzf->out, bgzf->outsize);
                break;
            }
            else if ((int)GetRCObject(rc) == rcTimeout
                     || (int)GetRCObject(rc) == rcData)
            {
                DBG("\t\tParser queue empty");
                if (KQueueSealed(que)) {
                    DBG("\t\tQueue sealed, Parser complete");
                    return false;
                }
            }
            else
            {
                ERR("Parser rc=%d", rc);
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            }
        }
        return true;
    }

  public:
    bool getbytes(KQueue* que, char* dest, size_t len)
    {
        DBG("Getting %d", len);
        if (len == 0) ERR("Empty get");
        while (len)
        {
            if (bgzf == NULL || bgzf->outsize == 0) {
                DBG("need %d more", len);
                if (!getnextBGZF(que)) return false;
                cur = bgzf->out;
            }

            size_t howmany = MIN(len, bgzf->outsize);

            memmove(dest, cur, howmany);
            len -= howmany;
            bgzf->outsize -= howmany;
            dest += howmany;
            cur += howmany;
        }
        return true;
    }

  private:
    BGZF* bgzf;
    Bytef* cur;
};

static BGZFview bview;

static rc_t seeker(const KThread* kt, void* in)
{
    Extractor* state = (Extractor*)in;
    pthread_t threadid = pthread_self();
    DBG("\tSeeker thread %p %lu started.", kt, threadid);

    state->file_pos = 0;
    while (true)
    {
        if (state->readbuf_sz < 28) {
            ERR("Small block:%d", state->readbuf_sz);
        }

        if (!memcmp(state->readbuf,
                    "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43"
                    "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00"
                    "\x00",
                    28))
        {
            DBG("complete EOF marker found");
            break;
        }

        z_stream strm;
        memset(&strm, 0, sizeof strm);
        strm.next_in = (Bytef*)state->readbuf;
        strm.avail_in = (uInt)state->readbuf_sz;
        int zrc = inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
        switch (zrc)
        {
        case Z_OK:
            break;
        case Z_MEM_ERROR:
            ERR("error: Out of memory in zlib");
            return RC(rcAlign, rcFile, rcReading, rcMemory, rcExhausted);
        case Z_VERSION_ERROR:
            ERR("zlib version is not compatible; need version %s but "
                "have %s",
                ZLIB_VERSION, zlibVersion());
            return RC(rcAlign, rcFile, rcConstructing, rcNoObj, rcUnexpected);
        case Z_STREAM_ERROR:
            ERR("zlib stream error");
            return RC(rcAlign, rcFile, rcConstructing, rcNoObj, rcUnexpected);
        default:
            ERR("zlib error");
            return RC(rcAlign, rcFile, rcConstructing, rcNoObj, rcUnexpected);
        }

        gz_header head;
        u8 extra[256];
        memset(&head, 0, sizeof head);
        head.extra = extra;
        head.extra_max = sizeof(extra);
        char outbuf[64];
        strm.next_out = (Bytef*)outbuf;
        strm.avail_out = 64;
        zrc = inflateGetHeader(&strm, &head);
        while (head.done == 0)
        {
            DBG("inflating gzip header");
            zrc = inflate(&strm, Z_BLOCK);
            if (zrc != Z_OK) {
                for (int i = 0; i != 4; ++i)
                    DBG("readbuf: %x", (unsigned char)state->readbuf[i]);
                ERR("inflate error %d %s", zrc, strm.msg);
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            }
        }

        DBG("found gzip header");
        // BC 02 bb
        if (head.extra && head.extra_len == 6 && head.extra[0] == 'B'
            && head.extra[1] == 'C' && head.extra[2] == 2
            && head.extra[3] == 0)
        {
            u16 bsize = head.extra[4] + head.extra[5] * 256;
            inflateEnd(&strm);
            DBG("total_in:%d", strm.avail_in);
            if (bsize <= 28) {
                DBG("small block found");
                break;
            }

            size_t block_size = 12; // Up to and including XLEN
            block_size += head.extra_len;
            block_size += (bsize - head.extra_len - 19); // CDATA
            block_size += 8;                             // CRC32 and isize
            DBG("block_size is %d bsize is %d", block_size, bsize);

            state->file_pos += block_size;

            BGZF* bgzf = (BGZF*)calloc(1, sizeof(BGZF));
            KLockMake(&bgzf->lock);
            KLockAcquire(bgzf->lock); // Not ready for parsing
            bgzf->state = compressed;
            bgzf->insize = bsize + 1;
            memmove(bgzf->in, state->readbuf, block_size + 1);
            bgzf->outsize = sizeof(bgzf->out);

            struct timeout_t tm;
            TimeoutInit(&tm, 1000); // 1 second
            while (true)
            {
                // Add to Inflate queue
                rc_t rc = KQueuePush(state->inflatequeue, (void*)bgzf, &tm);
                if ((int)GetRCObject(rc) == rcTimeout) {
                    DBG("inflate queue full");
                }
                else if (rc == 0)
                {
                    DBG("inflate queued: %p %d %d %d", bgzf->in, bgzf->insize,
                        rc, rcTimeout);
                    break;
                }
                else
                {
                    ERR("inflate queue %d", rc);
                }
            }

            while (true)
            {
                // Add to parse queue
                // lock will prevent parsing until inflater
                // thread finished with this chunk.
                rc_t rc = KQueuePush(state->parsequeue, (void*)bgzf, &tm);
                if ((int)GetRCObject(rc) == rcTimeout) {
                    DBG("parse queue full");
                }
                else if (rc == 0)
                {
                    DBG("parse queued: %p %d %d %d", bgzf->in, bgzf->insize,
                        rc, rcTimeout);
                    break;
                }
                else
                {
                    DBG("parse queued%d", rc);
                }
            }
        }
        else
        {
            ERR("error: BAM required extra extension not found");
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }

        DBG("reading in at %d", state->file_pos);
        rc_t rc = KFileReadAll(state->infile, state->file_pos, state->readbuf,
                               state->readbuf_sz, &state->readbuf_sz);
        // state->file_pos+=state->readbuf_sz;
        if (rc) {
            ERR("readfile error");
            return rc;
        }
        DBG("Read in %d", state->readbuf_sz);
        state->readbuf_pos = 0;
        //        if (state->file_pos > 10000000) state->readbuf_sz = 0;
        //        TODO
        if (state->readbuf_sz == 0) {
            DBG("Buffer complete. EOF");
            break;
        }
    }

    DBG("seeker thread complete");
    KQueueSeal(state->parsequeue);
    return 0;
}

static rc_t inflater(const KThread* kt, void* in)
{
    Extractor* state = (Extractor*)in;
    struct timeout_t tm;

    z_stream strm;
    pthread_t threadid = pthread_self();
    DBG("\tInflater thread %p %lu started.", kt, threadid);

    while (true)
    {
        void* where = NULL;
        DBG("\t\tthread %lu checking queue", threadid);
        TimeoutInit(&tm, 1000); // 1 seconds
        rc_t rc = KQueuePop(state->inflatequeue, &where, &tm);
        if (rc == 0) {
            BGZF* bgzf = (BGZF*)where;
            DBG("\t\tinflater thread %lu BGZF %p size %u", threadid, bgzf->in,
                bgzf->insize);
            if (bgzf->state != compressed) {
                ERR("Inflater bad state");
                return RC(rcAlign, rcFile, rcReading, rcData, rcInvalid);
            }

            memset(&strm, 0, sizeof strm);
            DBG("\tinflating %d bytes", bgzf->insize);
            strm.next_in = bgzf->in;
            strm.avail_in = bgzf->insize;
            strm.next_out = bgzf->out;
            strm.avail_out = bgzf->outsize;
            int zrc = inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
            switch (zrc)
            {
            case Z_OK:
                break;
            case Z_MEM_ERROR:
                ERR("Out of memory in zlib");
                return RC(rcAlign, rcFile, rcReading, rcMemory, rcExhausted);
            case Z_VERSION_ERROR:
                ERR("zlib version is not compatible; need "
                    "version %s but "
                    "have %s",
                    ZLIB_VERSION, zlibVersion());
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            case Z_STREAM_ERROR:
                ERR("zlib stream error");
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            default:
                ERR("zlib error %s", strm.msg);
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            }

            zrc = inflate(&strm, Z_FINISH);
            switch (zrc)
            {
            case Z_OK:
            case Z_STREAM_END:
                DBG("\t\tthread %lu OK %d %d %lu", threadid, strm.avail_in,
                    strm.avail_out, strm.total_out);
                bgzf->outsize = strm.total_out;
                bgzf->state = uncompressed;
                DBG("Ready for parsing, unlocking");
                KLockUnlock(bgzf->lock); // OK to parse now
                break;
            case Z_MEM_ERROR:
                ERR("error: Out of memory in zlib");
                return RC(rcAlign, rcFile, rcReading, rcMemory, rcExhausted);
            case Z_VERSION_ERROR:
                ERR("zlib version is not compatible; need "
                    "version %s but "
                    "have %s",
                    ZLIB_VERSION, zlibVersion());
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            case Z_STREAM_ERROR:
                ERR("zlib stream error %s", strm.msg);
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            default:
                ERR("inflate error %d %s", zrc, strm.msg);
                return RC(rcAlign, rcFile, rcConstructing, rcNoObj,
                          rcUnexpected);
            }
            inflateEnd(&strm);
        }
        else if ((int)GetRCObject(rc) == rcTimeout
                 || (int)GetRCObject(rc) == rcData)
        {
            DBG("\t\tthread %lu queue empty", threadid);
            if (KQueueSealed(state->parsequeue)) {
                DBG("\t\tqueue sealed, inflater thread %lu terminating.",
                    threadid);
                return 0;
            }
        }
        else
        {
            ERR("rc=%d", rc);
            return rc;
        }
    }

    ERR("\t\tinflater thread %lu wrongly terminating.", threadid);
    return 0;
}

rc_t BAMGetHeaders(Extractor* state)
{
    DBG("BAMGetHeaders");
    char magic[4];
    if (!bview.getbytes(state->parsequeue, magic, 4)) return 1;
    if (memcmp(magic, "BAM\x01", 4)) {
        ERR("BAM magic not found");
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }
    i32 l_text;
    if (!bview.getbytes(state->parsequeue, (char*)&l_text, 4)) return 1;
    if (l_text < 0) {
        ERR("error: invalid l_text");
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }

    char* text = (char*)pool_calloc(l_text + 2);
    if (!bview.getbytes(state->parsequeue, text, l_text)) return 1;

    DBG("SAM header %d %d:'%s'", l_text, strlen(text), text);
    char* t = text;
    while (strlen(t))
    {
        char* nl = (char*)strchr(t, '\n');
        if (!nl) {
            size_t linelen = strlen(t);
            DBG("noln linelen %d", linelen);
            memmove(curline, t, linelen);
            curline[linelen + 1] = '\n';
            curline[linelen + 2] = '\0';
            t += linelen;
        }
        else
        {
            size_t linelen = 1 + nl - t;
            DBG("ln   linelen %d", linelen);
            memmove(curline, t, linelen);
            curline[linelen] = '\0';
            t += linelen;
        }
        DBG("curline is '%s'", curline);
        if (curline[0] != '@') ERR("Not a SAM header line: '%s'", curline);
        curline_len = strlen(curline);
        SAMparse(state);
    }
    pool_free(text);
    text = NULL;

    if (!bview.getbytes(state->parsequeue, (char*)&state->n_ref, 4)) return 1;
    if (state->n_ref < 0) {
        ERR("error: invalid n_ref: %d", state->n_ref);
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }
    DBG("# references %d", state->n_ref);

    for (int i = 0; i != state->n_ref; ++i)
    {
        i32 l_name;
        if (!bview.getbytes(state->parsequeue, (char*)&l_name, 4)) return 1;
        DBG("ref #%d/%d: l_name=%d", i, state->n_ref, l_name);
        if (l_name < 0) {
            ERR("error: invalid reference name length:%d", l_name);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }
        if (l_name > 256) {
            ERR("warning: Long reference name:%d", l_name);
            //            return 0;
        }
        char* name = (char*)pool_alloc(l_name + 1);
        if (!bview.getbytes(state->parsequeue, name, l_name)) return 1;
        DBG("ref #%d:'%s'", i, name);
        pool_free(name);
        name = NULL;
        i32 l_ref;
        if (!bview.getbytes(state->parsequeue, (char*)&l_ref, 4)) return 1;
        DBG("length of reference sequence %d=%d", i, l_ref);
    }
    DBG("End of references");

    DBG("BAM done with headers");

    return 0;
}

rc_t BAMGetAlignments(Extractor* state)
{
    bamalign align;

    while (bview.getbytes(state->parsequeue, (char*)&align, sizeof(align)))
    {
        DBG("alignment block_size=%d refID=%d pos=%d", align.block_size,
            align.refID, align.pos);

        if (align.block_size < 0) {
            ERR("error: invalid block_size:%d", align.block_size);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }
        if (align.pos < -1) {
            ERR("error: invalid pos:%d", align.pos);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }
        if (align.refID < -1 || align.refID > state->n_ref) {
            ERR("error: bad refID:%d", align.refID);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }
        if (align.next_refID < -1 || align.next_refID > state->n_ref) {
            ERR("error: bad next_refID:%d", align.next_refID);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }
        if (align.next_pos < -1) {
            ERR("error: bad next_pos:%d", align.next_pos);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }

        if (align.tlen < 0) {
            ERR("error: bad tlen:%d", align.tlen);
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        }

        DBG("align.bin_mq_nl=%d", align.bin_mq_nl);
        u16 bin = align.bin_mq_nl >> 16;
        u8 mapq = (align.bin_mq_nl >> 8) & 0xff;
        u8 l_read_name = align.bin_mq_nl & 0xff;
        DBG("bin=%d mapq=%d l_read_name=%d", bin, mapq, l_read_name);
        if (l_read_name > 64) ERR("Long (%d) read_name", l_read_name);

        u16 flag = align.flag_nc >> 16;
        u16 n_cigar_op = align.flag_nc & 0xffff;
        DBG("flag=%x n_cigar_op=%d", flag, n_cigar_op);

        char* read_name = (char*)pool_alloc(l_read_name);
        if (!bview.getbytes(state->parsequeue, read_name, l_read_name))
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        DBG("read_name='%s'", read_name);

        char* scigar = NULL;
        if (n_cigar_op > 0) {
            static const char opmap[] = "MIDNSHP=X???????";
            u32* cigar = (u32*)pool_alloc(n_cigar_op * sizeof(u32));

            if (!bview.getbytes(state->parsequeue, (char*)cigar,
                                n_cigar_op * 4))
                return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
            // Compute size of expanded CIGAR
            size_t rleopslen = 0;
            for (int i = 0; i != n_cigar_op; ++i)
            {
                i32 oplen = cigar[i] >> 4;
                if (!oplen) ERR("Bogus CIGAR op length");
                rleopslen += oplen;
            }
            DBG("rleopslen %d", rleopslen);
            scigar = (char*)pool_calloc(rleopslen + 1);
            char* p = scigar;
            for (int i = 0; i != n_cigar_op; ++i)
            {
                i32 oplen = cigar[i] >> 4;
                i32 op = cigar[i] & 0xf;
                DBG("\tcigar %d=%x len=%d %d(%c)", i, cigar[i], oplen, op,
                    opmap[op]);
                for (int j = 0; j != oplen; ++j)
                    *p++ = (char)opmap[op];
            }
            *p = '\0';
            pool_free(cigar);
            cigar = NULL;
        }
        else
        {
            scigar = pool_strdup("");
        }
        DBG("scigar is '%s'", scigar);

        static const char seqmap[] = "=ACMGRSVTWYHKDBN";
        char* seq = (char*)pool_calloc(align.l_seq + 1);
        int bytesofseq = (align.l_seq + 1) / 2;
        char* seqbytes = (char*)pool_alloc(bytesofseq);
        if (!bview.getbytes(state->parsequeue, seqbytes, bytesofseq))
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
        int i = 0;
        int j = 0;
        while (i < align.l_seq)
        {
            seq[i] = seqmap[seqbytes[j] >> 4];
            i += 2;
            j += 1;
        }
        i = 1;
        j = 0;
        while (i < align.l_seq)
        {
            seq[i] = seqmap[seqbytes[j] & 0xf];
            i += 2;
            j += 1;
        }

        char* qual = (char*)pool_alloc(align.l_seq);
        if (!bview.getbytes(state->parsequeue, qual, align.l_seq))
            return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);

        DBG("%d pairs in sequence", align.l_seq);
        DBG("seq='%s'", seq);

        int remain = align.block_size
                     - (sizeof(align) + l_read_name + n_cigar_op * 4
                        + bytesofseq + align.l_seq) + 4; // TODO, why 4?
        DBG("%d bytes remaining for ttvs", remain);
        char* ttvs = NULL;
        if (remain) {
            ttvs = (char*)pool_alloc(remain);
            if (!bview.getbytes(state->parsequeue, ttvs, remain))
                return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
            char* cur = ttvs;
            while (cur < ttvs + remain)
            {
                char tag[2];
                char c;
                i8 i8;
                u8 u8;
                i16 i16;
                u16 u16;
                i32 i32;
                u32 u32;
                char* z;
                tag[0] = *cur++;
                tag[1] = *cur++;
                char val_type = *cur++;
                DBG("ttv: %c%c:%c", tag[0], tag[1], val_type);
                switch (val_type)
                {
                case 'A':
                    c = *cur++;
                    DBG("val='%c'", c);
                    break;
                case 'c':
                    i8 = *cur++;
                    DBG("val=%d", i8);
                    break;
                case 'C':
                    u8 = *cur++;
                    DBG("val=%d", u8);
                    break;
                case 's':
                    memmove(&i16, cur, 2);
                    DBG("val=%d", i16);
                    cur += 2;
                    break;
                case 'S':
                    memmove(&u16, cur, 2);
                    DBG("val=%d", u16);
                    cur += 2;
                    break;
                case 'i':
                    memmove(&i32, cur, 4);
                    cur += 4;
                    break;
                case 'I':
                    memmove(&u32, cur, 4);
                    cur += 4;
                    break;
                case 'f':
                    // float f;
                    break;
                case 'Z':
                    z = cur;
                    while (isprint(*cur))
                        ++cur;
                    DBG("val='%s'", z);
                    ++cur;
                    break;
                case 'H':
                    z = cur;
                    while (isalnum(*cur))
                        ++cur;
                    DBG("val='%s'", z);
                    // TODO: Convert to ?
                    ++cur;
                    break;
                case 'B':
                    val_type = *cur++;
                    memmove(&u32, cur, 4);
                    cur += 4;
                    cur += u32 * 1; // TODO, based on size of
                                    // val_type
                    break;
                default:
                    ERR("Bad val_type:%c", val_type);
                    return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
                }
            }
        }
        DBG("no more ttvs");

        // We want read (seq), cigar, rname, pos and flags
        // name=Qname, reference sequence name
        // read_name=qname, query template name
        process_alignment(state,
                          NULL, // QNAME
                          flag,
                          read_name, // RNAME
                          align.pos,
                          "0",    // mapq
                          scigar, // cigar
                          NULL,   // rnext
                          "0",    // pnext
                          "0",    // tlen
                          seq,    // read
                          NULL    // qual
                          );

        pool_free(read_name);
        read_name = NULL;
        pool_free(seq);
        seq = NULL;
        pool_free(seqbytes);
        seqbytes = NULL;
        pool_free(qual);
        qual = NULL;
        pool_free(ttvs);
        ttvs = NULL;
        pool_free(scigar);
        scigar = NULL;

        if (VectorLength(&state->alignments) == 64) {
            DBG("Have %d BAM alignments", VectorLength(&state->alignments));
            return 0;
        }
    }

    if (!KQueueSealed(state->parsequeue)) {
        ERR("out of data but queue not sealed");
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }
    DBG("parser complete");

    return 0;
}

rc_t releasethreads(Extractor* state)
{
    rc_t rc;
    DBG("Releasing threads");
    for (u32 i = 0; i != VectorLength(&state->threads); ++i)
    {
        KThread* t = (KThread*)VectorGet(&state->threads, i);
        rc = KThreadRelease(t);
        if (rc) return rc;
    }
    DBG("Released threads");
    return 0;
}

rc_t threadinflate(Extractor* state)
{
    rc_t rc;
    DBG("Starting threads");

    // Inflater threads
    for (i32 i = 0; i != MAX(1, state->num_threads - 1); ++i)
    {
        KThread* inflaterthread;
        rc = KThreadMake(&inflaterthread, inflater, (void*)state);
        if (rc) return rc;
        rc = KThreadDetach(inflaterthread);
        if (rc) return rc;
        VectorAppend(&state->threads, NULL, inflaterthread);
    }

    // Seeker thread
    KThread* seekerthread;
    rc = KThreadMake(&seekerthread, seeker, (void*)state);
    if (rc) return rc;
    rc = KThreadDetach(seekerthread);
    if (rc) return rc;
    VectorAppend(&state->threads, NULL, seekerthread);

    DBG("Threads started.");

    return 0;
}
