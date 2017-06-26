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

#include "/home/vartanianmh/devel/sra-tools/tools/general-loader/general-writer.hpp"

#include <align/samextract-lib.h>
#include <ctype.h>
#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/directory.h>
#include <kfs/file.h>
#include <klib/defs.h>
#include <klib/rc.h>
#include <klib/text.h>
#include <klib/vector.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

extern "C" {
ver_t CC KAppVersion(void) { return 0x1000000; }

rc_t CC UsageSummary(char const* name)
{
    fprintf(stderr, "Usage: %s file.{sb}am [file...]\n", name);
    return 0;
}

rc_t CC Usage(Args const* args) { return 0; }

rc_t CC KMain(int argc, char* argv[])
{
    rc_t rc;
    if (argc == 1) {
        UsageSummary(argv[0]);
        return 0;
    }
    while (--argc) {
        const char* fname = *(++argv);

        struct KDirectory* srcdir = NULL;
        const struct KFile* infile = NULL;
        rc = KDirectoryNativeDir(&srcdir);
        if (rc) return rc;

        rc = KDirectoryOpenFileRead(srcdir, &infile, fname);
        KDirectoryRelease(srcdir);
        if (rc) return rc;
        srcdir = NULL;

        String sfname;
        StringInitCString(&sfname, fname);
        struct timespec stime, etime;
        clock_gettime(CLOCK_REALTIME, &stime);
        SAMExtractor* extractor;
        rc = SAMExtractorMake(&extractor, infile, &sfname, -1);
        fprintf(stderr, "Made extractor for %s\n", fname);
        if (rc) return rc;

        //        rc = SAMExtractorAddFilter(extractor, NULL, -1, -1, false);
        if (rc) return rc;

        Vector headers;
        rc = SAMExtractorGetHeaders(extractor, &headers);
        if (rc) return rc;
        fprintf(stderr, "Got %d headers\n", VectorLength(&headers));
#if 0
        GeneralWriter * gw=new GeneralWriter(1, 32768);
        GeneralWrite &out=*gw;
        gw->setSoftwarename("samextract 0.1");
        gw->setRemotePath("samextract.db");
        gw->useSchema("bamdb.schema");
        tbl_id=gw->addTable("header");
        keyid=gw->addcolumn("hdrkey");
#endif

        for (uint32_t i = 0; i != VectorLength(&headers); ++i) {
            Header* hdr = (Header*)VectorGet(&headers, i);
            Vector* tvs = &hdr->tagvalues;
            //            fprintf(stderr,"\tHeader%d: %s\n", i,
            //            hdr->headercode);
            for (uint32_t j = 0; j != VectorLength(tvs); ++j) {
                TagValue* tv = (TagValue*)VectorGet(tvs, j);

                //                fprintf(stderr,"\t\t%d\t%s %s\n", j,
                //                tv->tag, tv->value);
            }
            // Do stuff with headers
        }
        SAMExtractorInvalidateHeaders(extractor);

        fprintf(stderr, "Getting Alignments\n");
        int total = 0;
        uint32_t vlen;
        do {
            Vector alignments;
            rc = SAMExtractorGetAlignments(extractor, &alignments);
            if (rc) {
                fprintf(stderr, "GetAligned returned rc\n");
                return rc;
            }
            vlen = VectorLength(&alignments);
            total += vlen;
            //            fprintf(stderr, "Got %d alignments\n", total);
            //            fprintf(stderr,"\n\nReturned %d alignments\n",vlen);
            for (uint32_t i = 0; i != vlen; ++i) {
                Alignment* align = (Alignment*)VectorGet(&alignments, i);
                //                if (strlen(align->cigar) > 0)
                //                    fprintf(stderr, "cigar is %s\n",
                //                    align->cigar);
                //                fprintf(stderr,"\tAlignment%2d: %s\n", i,
                //                align->read);
                // Do stuff with headers
            }
            //            fprintf(stderr,"\n");
            SAMExtractorInvalidateAlignments(extractor);
            // if (total > 100000) break;
        } while (vlen);

        rc = SAMExtractorRelease(extractor);
        if (rc) {
            fprintf(stderr, "ExtractorRelease returned rc %d\n", rc);
            return rc;
        }

        fprintf(stderr, "Done with file, %d alignments\n", total);
        clock_gettime(CLOCK_REALTIME, &etime);
        u64 nanos = etime.tv_sec - stime.tv_sec;
        nanos *= 1000000000;
        nanos += (etime.tv_nsec - stime.tv_nsec);
        fprintf(stderr, "Parse time %lu ms", nanos / 1000000l);

        KFileRelease(infile);
        infile = NULL;
    }
    return 0;
}
}
