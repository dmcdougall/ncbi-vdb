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

#ifndef _h_samextract_lib_
#define _h_samextract_lib_
#include <align/extern.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <klib/text.h>
#include <kproc/queue.h>
#include <kfs/file.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef enum file_type
{
    unknown,
    SAM,
    BAM,
    SAMGZUNSUPPORTED
} file_type;
//    typedef enum file_status { init, headers, alignments, done} file_status;

typedef struct SAMExtractor
{
    const KFile* infile;
    String* fname;

    Vector headers;
    Vector alignments;

    Vector tagvalues;
    Vector* prev_headers;
    Vector* prev_aligns;

    //        Vector chunks;
    int32_t num_threads;
    Vector threads;
    KQueue* inflatequeue;
    KQueue* parsequeue;

    uint32_t pos;
    uint64_t file_pos;

    char* readbuf;
    size_t readbuf_sz;
    size_t readbuf_pos;

    file_type file_type;

    int32_t n_ref;

    rc_t rc;

    //        file_status file_status;
    bool hashdvn;
    bool hashdso;
    bool hashdgo;
    bool hassqsn;
    bool hassqln;
    bool hasrgid;
    bool haspgid;
} SAMExtractor;

typedef struct tagvalue
{
    const char* tag; /* VN, SN, LN, ID, ... */
    const char* value;
} TagValue;

typedef struct Header
{
    const char* headercode; /* HD, SQ, RG, PG, CO */
    Vector tagvalues;
} Header;

typedef struct Alignment
{
    const char* read;
    const char* cigar;
    const char* rname;
    uint32_t pos;
    uint16_t flags;
} Alignment;

/* TODO: API change: Pass in filename for diagnostics */
ALIGN_EXTERN rc_t CC SAMExtractorMake(SAMExtractor** state, const KFile* fin,
                                      String* fname_desc,
                                      int32_t num_threads);
ALIGN_EXTERN rc_t CC SAMExtractorRelease(SAMExtractor* state); /* dtor */

ALIGN_EXTERN rc_t CC
    SAMExtractorGetHeaders(SAMExtractor* state, Vector* headers);
ALIGN_EXTERN rc_t CC SAMExtractorInvalidateHeaders(SAMExtractor* state);

ALIGN_EXTERN rc_t CC
    SAMExtractorGetAlignments(SAMExtractor* state, Vector* alignments);
ALIGN_EXTERN rc_t CC SAMExtractorInvalidateAlignments(SAMExtractor* state);

#ifdef __cplusplus
}
#endif
#endif /* __h_sam_extract_lib_ */
