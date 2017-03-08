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

#ifndef _h_samextract_lib_
#define _h_samextract_lib_
#include <align/extern.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct Extractor
    {
        char *mmapbuf;
        off_t mmapbuf_sz;
        char *mmapbuf_cur;

        Vector headers;
        Vector alignments;
        Vector tagvalues; // temp

        Vector allocs;

        Vector * prev_headers;
        Vector * prev_aligns;

        char * read;
        char * cigar;
        char * rname;
        uint32_t pos;

        char * tags; // Space delimited tags seen in current line
        char * seqnames;
        char * ids;
        rc_t rc;
    } Extractor;

    typedef struct tagvalue
    {
        const char * tag; // VN, SN, LN, ID, ...
        const char * value;
    } TagValue;

    typedef struct Header
    {
        const char * headercode; // HD, SQ, RG, PG, CO
        Vector tagvalues;
    } Header;

    typedef struct Alignment
    {
        const char * read;
        const char * cigar;
        const char * rname;
        uint32_t pos;
    } Alignment;

    ALIGN_EXTERN rc_t CC SAMExtractorMake(Extractor **state, const char * fname, uint32_t num_threads);
    ALIGN_EXTERN rc_t CC SAMExtractorRelease(Extractor *state); // dtor

    ALIGN_EXTERN rc_t CC SAMExtractorGetHeaders(Extractor *state, Vector *headers);
    ALIGN_EXTERN rc_t CC SAMExtractorInvalidateHeaders(Extractor *state);

    ALIGN_EXTERN rc_t CC SAMExtractorGetAlignments(Extractor *state, Vector *alignments);
    ALIGN_EXTERN rc_t CC SAMExtractorInvalidateAlignments(Extractor *state);

#ifdef __cplusplus
}
#endif
#endif // __h_sam_extract_lib_

