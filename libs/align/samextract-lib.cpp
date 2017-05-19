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

#define __STDC_LIMIT_MACROS
#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/file.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <kproc/queue.h>
#include <kproc/thread.hpp>
#include <kproc/timeout.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <strtol.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "samextract.h"
#include "samextract-tokens.h"
#include "samextract-pool.h"
#include <align/samextract-lib.h>
#include <stdint.h>

char curline[READBUF_SZ + 1];
size_t curline_len = 0;
String* fname_desc = NULL;

void logmsg(const char* fname, int line, const char* func,
            const char* severity, const char* fmt, ...)
{
    char* buf;
    size_t bufsize = 0;
    FILE* buffd;

    pthread_t threadid = pthread_self();

    const char* basename = strrchr(fname, '/');
    if (!basename) basename = strrchr(fname, '\\');
    if (basename) ++basename;
    if (!basename) basename = fname;
    va_list args;
    va_start(args, fmt);

    buffd = open_memstream(&buf, &bufsize);
    if (buffd == NULL) {
        fprintf(stderr, "can't open memstream\n");
        return;
    }
    fprintf(buffd, "%s(%lu) ", severity, threadid % 100);
    if (fname_desc) fprintf(buffd, "`%s`:", fname_desc->addr);
    vfprintf(buffd, fmt, args);
    va_end(args);
    fprintf(buffd, "\t[%s:%s():%d]\n", basename, func, line);
    fclose(buffd);
    size_t r = fwrite(buf, bufsize, 1, stderr);
    if (r != 1)
        fprintf(stderr, "previous %zd log message truncated\n", bufsize);
    free(buf);
    buf = NULL;
    fflush(stderr);
    if (!strcmp(severity, "Error")) abort();
}

rc_t SAM_parseline(SAMExtractor* state)
{
    state->rc = 0;
    DBG("Parsing line (%d bytes): '%s'", strlen(curline), curline);
    SAMparse(state);
    return state->rc;
}

int moredata(char* buf, int* numbytes, size_t maxbytes)
{
    if (!curline_len)
        DBG("nomoredata");
    else
        DBG("  moredata %p %d\ncurline:'%s'", buf, maxbytes, curline);
    memmove(buf, curline, curline_len);
    *numbytes = curline_len;
    curline_len = 0;
    return 0;
}

inline rc_t readfile(SAMExtractor* state)
{
    if (state->readbuf_pos == state->readbuf_sz) {
        state->readbuf_sz = READBUF_SZ;
        DBG("reading in at %d", state->file_pos);
        rc_t rc = KFileReadAll(state->infile, state->file_pos, state->readbuf,
                               state->readbuf_sz, &state->readbuf_sz);
        state->file_pos += state->readbuf_sz;
        if (rc) {
            ERR("readfile error");
            return rc;
        }
        DBG("Read in %d", state->readbuf_sz);
        state->readbuf_pos = 0;
        if (!state->readbuf_sz) {
            DBG("Buffer complete. EOF");
        }
    }
    return 0;
}

void SAMerror(SAMExtractor* state, const char* s)
{
    ERR(" Parsing error: %s\nLine was:'%s'", s, curline);
    rc_t rc = RC(rcAlign, rcRow, rcParsing, rcData, rcInvalid);
    state->rc = rc;
    abort();
}

/* low<=str<=high */
bool inrange(const char* str, i64 low, i64 high)
{
    i64 i = strtoi64(str, NULL, 10);
    if (errno) return false;
    if (i < low || i > high) return false;
    return true;
}

bool ismd5(const char* str)
{
    size_t i;
    size_t len = strlen(str);

    if (len != 32) return false;

    for (i = 0; i != len; ++i)
    {
        if (!isalnum(str[i]) && str[i] != '*') return false;
    }

    return true;
}

/* Avoiding handling this as flex token because so many other ReadGroup values
 * could end up looking like flow orders.
 */
bool isfloworder(const char* str)
{
    size_t i;
    size_t len = strlen(str);

    if (len == 1 && str[0] == '*') return true;
    for (i = 0; i != len; ++i)
    {
        switch (str[i])
        {
        case 'A':
        case 'C':
        case 'M':
        case 'G':
        case 'R':
        case 'S':
        case 'V':
        case 'T':
        case 'W':
        case 'Y':
        case 'H':
        case 'K':
        case 'D':
        case 'B':
        case 'N':
            continue;
        default:
            return false;
        }
    }
    return true;
}

rc_t process_header(SAMExtractor* state, const char* type, const char* tag,
                    const char* value)
{
    DBG("processing type:%s tag:%s value:%s", type, tag, value);
    if (strcmp(type, "HD") && strcmp(type, "SQ") && strcmp(type, "RG")
        && strcmp(type, "PG"))
    {
        ERR("record '%s' must be HD, SQ, RG or PG", type);
        rc_t rc = RC(rcAlign, rcRow, rcParsing, rcData, rcInvalid);
        state->rc = rc;
        return rc;
    }

    if (strlen(tag) != 2) {
        ERR("tag '%s' must be 2 characters", tag);
        rc_t rc = RC(rcAlign, rcRow, rcParsing, rcData, rcInvalid);
        state->rc = rc;
        return rc;
    }

    if (islower(tag[0] && islower(tag[1]))) {
        DBG("optional tag");
    }

    TagValue* tv = (TagValue*)pool_alloc(sizeof(TagValue));
    tv->tag = pool_strdup(tag);
    tv->value = pool_strdup(value);
    VectorAppend(&state->tagvalues, NULL, tv);

    return 0;
}

rc_t mark_headers(SAMExtractor* state, const char* type)
{
    DBG("mark_headers");
    Header* hdr = (Header*)pool_alloc(sizeof(Header));
    hdr->headercode = type;
    VectorCopy(&state->tagvalues, &hdr->tagvalues);
    VectorAppend(&state->headers, NULL, hdr);
    VectorWhack(&state->tagvalues, NULL, NULL);
    return 0;
}

// Returns true if we can skip record
bool filter(const SAMExtractor* state, const char* rname, ssize_t pos)
{
    if (state->filter_rname && strcmp(state->filter_rname, rname))
        return true;

    if (state->file_type == SAM) --pos; // Internally use and expose 0-based

    if (pos < 0)
        return false; // No filtering if pos uncertain (0 in SAM, -1 in BAM)

    if (state->filter_pos != -1) {
        if (pos < state->filter_pos) // Before pos
            return true;

        if (state->filter_length != -1)
            if (pos > (state->filter_pos
                       + state->filter_length)) // After pos+length
                return true;
    }

    return false;
}

rc_t process_alignment(SAMExtractor* state, const char* qname, u16 flag,
                       const char* rname, i32 pos, const char* mapq,
                       const char* cigar, const char* rnext,
                       const char* pnext, const char* tlen, const char* seq,
                       const char* qual)
{
    DBG("process_alignment %s %d %u %s %s %s", rname, pos, flag, qname, rnext,
        qual);

    if (pos < -1) ERR("POS not in range %d", pos);

    if (state->file_type == SAM) {
        if (!inrange(mapq, 0, UINT8_MAX)) ERR("MAPQ not in range %s", mapq);

        if (!inrange(pnext, 0, INT32_MAX))
            ERR("PNEXT not in range %s", pnext);

        if (!inrange(tlen, INT32_MIN, INT32_MAX))
            ERR("TLEN not in range %s", tlen);
    }

    // TODO: cigar/rleopslen should be equal to l_seq

    // TODO: ordered
    if (filter(state, rname, pos)) {
        DBG("Skipping");
        return 0;
    }

    Alignment* align = (Alignment*)pool_alloc(sizeof(Alignment));

    align->read = seq;
    align->cigar = cigar;
    align->rname = rname;
    align->pos = pos;
    align->flags = flag;
    VectorAppend(&state->alignments, NULL, align);

    return 0;
}

// Reads next line into curline, returns false if file complete.
static bool readline(SAMExtractor* state)
{
    DBG("readline");
    if (readfile(state)) return false;
    char* line = curline;
    line[0] = '\0';
    curline_len = 0;
    // Is there a newline in current buffer?
    char* nl = (char*)memchr((state->readbuf + state->readbuf_pos), '\n',
                             (state->readbuf_sz - state->readbuf_pos));
    if (nl) {
        nl += 1;
        size_t len = nl - (state->readbuf + state->readbuf_pos);
        memmove(line, state->readbuf + state->readbuf_pos, len);
        curline_len += len;
        state->readbuf_pos += len;
        line[curline_len + 1] = '\0';
        return true;
    }

    // Nope, append and get more
    size_t len = (state->readbuf_sz - state->readbuf_pos);
    DBG("readline more %d/%d", state->readbuf_pos, len);
    memmove(line, state->readbuf + state->readbuf_pos, len);
    DBG("moreline was %d '%s'", strlen(line), line);
    line += len;
    curline_len += len;

    state->readbuf_pos = state->readbuf_sz;
    if (readfile(state)) return false;

    // Better be a newline now
    nl = (char*)memchr(state->readbuf, '\n', state->readbuf_sz);
    if (!nl) {
        return false;
    }
    DBG("found newline at %d", nl - state->readbuf);
    nl += 1;
    len = (nl - state->readbuf);
    memmove(line, state->readbuf, len);
    curline_len += len;
    state->readbuf_pos += len;
    line[curline_len + 1] = '\0';
    DBG("moreline  is %d %d '%s'", curline_len, strlen(line), line);

    return true;
}

LIB_EXPORT rc_t CC SAMExtractorMake(SAMExtractor** state, const KFile* fin,
                                    String* fname, int32_t num_threads = -1)
{
    SAMExtractor* s = (SAMExtractor*)malloc(sizeof(*s));
    *state = s;

    pool_init();

    s->infile = fin;
    s->fname = fname;
    fname_desc = fname;

    VectorInit(&s->headers, 0, 0);
    VectorInit(&s->alignments, 0, 0);
    VectorInit(&s->tagvalues, 0, 0);

    s->prev_headers = NULL;
    s->prev_aligns = NULL;

    s->num_threads = num_threads;

    // Default number of threads to number of cores
    if (s->num_threads <= -1)
        s->num_threads = sysconf(_SC_NPROCESSORS_ONLN) - 1;

    DBG("%d threads", s->num_threads);

    VectorInit(&s->threads, 0, 0);
    KQueueMake(&s->inflatequeue, 64);
    KQueueMake(&s->parsequeue, 64);

    s->pos = 0;
    s->file_pos = 0;

    s->readbuf = (char*)malloc(READBUF_SZ + 1);
    s->readbuf_sz = 0;
    s->readbuf_pos = 0;

    s->file_type = unknown;
    s->n_ref = -1;

    s->rc = 0;

    s->filter_rname = NULL;
    s->filter_pos = -1;
    s->filter_length = -1;
    s->filter_ordered = false;

    s->hashdvn = false;
    s->hashdso = false;
    s->hashdgo = false;
    s->hassqsn = false;
    s->hassqln = false;
    s->hasrgid = false;
    s->haspgid = false;

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorRelease(SAMExtractor* s)
{
    DBG("complete release_Extractor");
    SAMlex_destroy();

    SAMExtractorInvalidateHeaders(s);
    SAMExtractorInvalidateAlignments(s);

    releasethreads(s);

    pool_release();
    VectorWhack(&s->headers, NULL, NULL);
    KQueueRelease(s->inflatequeue);
    KQueueRelease(s->parsequeue);
    VectorWhack(&s->threads, NULL, NULL);
    free(s->readbuf);
    free(s->filter_rname);
    memset(s, 0, sizeof(SAMExtractor));
    free(s);

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorAddFilter(SAMExtractor* state,
                                         const char* rname, ssize_t pos,
                                         ssize_t length, bool ordered)
{
    // TODO: Check if GetHeaders/GetAlignments already invoked

    if (rname) state->filter_rname = strdup(rname);
    state->filter_pos = pos;
    state->filter_length = length;
    state->filter_ordered = ordered;

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorGetHeaders(SAMExtractor* s, Vector* headers)
{
    rc_t rc = 0;
    DBG("GetHeaders");

    u64 sz = 0;
    rc = KFileSize(s->infile, &sz);
    DBG("File size=%u", sz);
    if (sz < 12) {
        ERR("File too small");
        return RC(rcAlign, rcRow, rcParsing, rcData, rcInvalid);
    }

    rc = readfile(s);
    if (rc) return rc;

    if (!memcmp(s->readbuf, "\x1f\x8b\x08", 3)) {
        DBG("gzip file, BAM or SAM.gz");
        s->file_type = BAM;
    }
    else if (s->readbuf[0] == '@')
    {
        DBG("SAM file");
        s->file_type = SAM;
    }
    else
    {
        ERR("Unkown magic, not a SAM file.");
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }

    switch (s->file_type)
    {
    case BAM:
        rc = threadinflate(s);
        if (rc) return rc;
        rc = BAMGetHeaders(s);
        if (rc) return rc;
        break;
    case SAM:
        while (readline(s))
        {
            rc = SAM_parseline(s);
            if (rc) return rc;

            if (curline[0] != '@') {
                // First line of alignments will be processed
                DBG("out of headers");
                break;
            }
        }
        break;
    default:
        ERR("Unknown file type");
        return RC(rcAlign, rcFile, rcParsing, rcData, rcInvalid);
    }

    DBG("Done parsing headers");
    VectorInit(headers, 0, 0);
    VectorCopy(&s->headers, headers);
    s->prev_headers = headers;
    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorInvalidateHeaders(SAMExtractor* s)
{
    DBG("invalidate_headers");
    for (u32 i = 0; i != VectorLength(&s->headers); ++i)
    {
        Header* hdr = (Header*)VectorGet(&s->headers, i);

        hdr->headercode = NULL;

        Vector* tvs = &hdr->tagvalues;
        for (u32 j = 0; j != VectorLength(tvs); ++j)
        {
            TagValue* tv = (TagValue*)VectorGet(tvs, j);
            tv->tag = NULL;
            tv->value = NULL;
        }
        VectorWhack(&hdr->tagvalues, NULL, NULL);
        hdr = NULL;
    }
    pool_release();
    pool_init();
    VectorWhack(&s->headers, NULL, NULL);
    VectorWhack(&s->tagvalues, NULL, NULL);
    VectorWhack(s->prev_headers, NULL, NULL);
    s->prev_headers = NULL;
    return 0;
}

LIB_EXPORT rc_t CC
    SAMExtractorGetAlignments(SAMExtractor* s, Vector* alignments)
{
    rc_t rc = 0;
    SAMExtractorInvalidateAlignments(s);
    VectorInit(&s->alignments, 0, 0);
    VectorInit(alignments, 0, 0);

    if (s->file_type == SAM) {
        while (VectorLength(&s->alignments) < 64)
        {
            if (!readline(s)) break;

            if (curline[0] == '@') {
                ERR("header restarted");
                break;
            }

            rc = SAM_parseline(s);
            if (rc) return rc;
        }

        DBG("Done parsing %d alignments", VectorLength(&s->alignments));
    }
    else if (s->file_type == BAM)
    {
        rc = BAMGetAlignments(s);
        DBG("complete parsing %d alignments", VectorLength(&s->alignments));
        if (rc) {
            ERR("BAMGetAlignmentes failed");
            return rc;
        }
    }
    else
    {
        ERR("Unknown file type");
    }

    VectorCopy(&s->alignments, alignments);
    VectorWhack(&s->alignments, NULL, NULL);
    s->prev_aligns = alignments;

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorInvalidateAlignments(SAMExtractor* s)
{
    size_t num = VectorLength(&s->alignments);
    DBG("invalidate_alignments %d", num);
    pool_release();
    pool_init();
    VectorWhack(&s->alignments, NULL, NULL);
    VectorWhack(s->prev_aligns, NULL, NULL);
    s->prev_aligns = NULL;

    return 0;
}
