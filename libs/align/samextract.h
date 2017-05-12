/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnologmsgy Information
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

#ifndef _h_samextract_
#define _h_samextract_
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <kfs/file.h>
#include <kfs/buffile.h>
#include <zlib.h>
#include <sys/types.h>
#include <align/samextract-lib.h>

#define READBUF_SZ 65536

typedef int8_t i8;
typedef uint8_t u8;
typedef int16_t i16;
typedef uint16_t u16;
typedef int32_t i32;
typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;

typedef struct bamalign
{
    i32 block_size;
    i32 refID;
    i32 pos;
    u32 bin_mq_nl;
    u32 flag_nc;
    i32 l_seq;
    i32 next_refID;
    i32 next_pos;
    i32 tlen;
} bamalign;

#ifdef __cplusplus
extern "C" {
#endif
int SAMparse(SAMExtractor* state);
void SAMerror(SAMExtractor* state, const char* TODOmsg);
int SAMlex_destroy(void);
void logmsg(const char* fname, int line, const char* func,
            const char* severity, const char* fmt, ...);
int moredata(char* buf, int* numbytes, size_t maxbytes);
rc_t process_header(SAMExtractor* state, const char* type, const char* tag,
                    const char* value);
rc_t process_alignment(SAMExtractor* state, const char* qname, u16 flag,
                       const char* rname, i32 pos, const char* mapq,
                       const char* cigar, const char* rnext,
                       const char* pnext, const char* tlen, const char* seq,
                       const char* qual);
rc_t mark_headers(SAMExtractor* state, const char* type);
bool inrange(const char* str, i64 low, i64 high);
bool ismd5(const char* str);
bool isfloworder(const char* str);

void pool_init(void);
void pool_release(void);
void pool_free(void*);
void* pool_alloc(size_t sz);
void* pool_calloc(size_t sz);
char* pool_strdup(const char* str);
char* pool_memdup(const char* str, size_t len);
rc_t threadinflate(SAMExtractor* state);
rc_t BAMGetHeaders(SAMExtractor* state);
rc_t BAMGetAlignments(SAMExtractor* state);
rc_t releasethreads(SAMExtractor* state);
rc_t readfile(SAMExtractor* state);
extern char curline[];
extern size_t curline_len;
#ifdef __cplusplus
}
#endif
void samload(char const path[]);

#ifndef DEBUG
#define DEBUG 0
#endif

#define ERR(...) logmsg(__FILE__, __LINE__, __func__, "Error", __VA_ARGS__)

#define WARN(...) logmsg(__FILE__, __LINE__, __func__, "Warning", __VA_ARGS__)
#define INFO(...) logmsg(__FILE__, __LINE__, __func__, "Info", __VA_ARGS__)
#define DBG(...)                                                             \
    do                                                                       \
    {                                                                        \
        if (DEBUG)                                                           \
            logmsg(__FILE__, __LINE__, __func__, "Debug", __VA_ARGS__);      \
    } while (0)

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#endif
