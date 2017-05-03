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

#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/file.h>
#include <kfs/directory.h>
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
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <regex.h>
#include <stdint.h>
#include <unistd.h>
#include <align/samextract-lib.h>

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
    while (--argc)
    {
        const char* fname = *(++argv);

        struct KDirectory* srcdir = NULL;
        const struct KFile* infile = NULL;
        rc = KDirectoryNativeDir(&srcdir);
        if (rc) return rc;

        rc = KDirectoryOpenFileRead(srcdir, &infile, fname);
        KDirectoryRelease(srcdir);
        if (rc) return rc;
        srcdir = NULL;

        Extractor* extractor;
        rc_t rc = SAMExtractorMake(&extractor, infile, -1);
        fprintf(stderr, "Made extractor for %s\n", fname);
        if (rc) return rc;

        Vector headers;
        rc = SAMExtractorGetHeaders(extractor, &headers);
        if (rc) return rc;
        fprintf(stderr, "\n\nGot %d headers\n", VectorLength(&headers));
        for (uint32_t i = 0; i != VectorLength(&headers); ++i)
        {
            Header* hdr = (Header*)VectorGet(&headers, i);
            Vector* tvs = &hdr->tagvalues;
            //            fprintf(stderr,"\tHeader%d: %s\n", i,
            //            hdr->headercode);
            for (uint32_t j = 0; j != VectorLength(tvs); ++j)
            {
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
        do
        {
            Vector alignments;
            rc = SAMExtractorGetAlignments(extractor, &alignments);
            if (rc) {
                fprintf(stderr, "GetAligned returned rc\n");
                return rc;
            }
            vlen = VectorLength(&alignments);
            total += vlen;
fprintf(stderr,"Got %d alignments\n",total);
            //            fprintf(stderr,"\n\nReturned %d alignments\n",vlen);
            for (uint32_t i = 0; i != vlen; ++i)
            {
                Alignment* align = (Alignment*)VectorGet(&alignments, i);

                //                fprintf(stderr,"\tAlignment%2d: %s\n", i,
                //                align->read);
                // Do stuff with headers
            }
            //            fprintf(stderr,"\n");
            SAMExtractorInvalidateAlignments(extractor);
            // if (total > 100000) break;
        } while (vlen);

        SAMExtractorRelease(extractor);
        fprintf(stderr, "Done with file, %d alignments\n", total);

        KFileRelease(infile);
        infile = NULL;
    }
    fprintf(stderr, "KMain done\n");
    return 0;
}
