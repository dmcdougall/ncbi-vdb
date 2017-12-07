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
#include <linux/limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <vdb/blob.h>
#include <vdb/cursor.h>
#include <vdb/database.h>
#include <vdb/manager.h>
#include <vdb/table.h>

static const char tool_name[] = "irvrfy";

namespace ncbi
{
static rc_t process(const char* dbname)
{
    rc_t rc;

    KDirectory* srcdir = NULL;

    rc = KDirectoryNativeDir(&srcdir);
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    VDBManager* mgr = NULL;
    rc = VDBManagerMakeUpdate(&mgr, NULL); // NULL=No working directory
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    VDatabase* db = NULL;
    rc = VDBManagerOpenDBUpdate(mgr, &db, NULL, dbname);
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    const VTable* tbl = NULL;
    rc = VDatabaseOpenTableRead(db, &tbl, "hdrs");
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    const VCursor* curs = NULL;
    rc = VTableCreateCursorRead(tbl, &curs);
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }
    uint32_t group_idx = 0; // HDR, TAG, VALUE
    uint32_t hdr_idx = 0;
    uint32_t tag_idx = 0;
    uint32_t value_idx = 0;
    rc = VCursorAddColumn(curs, &group_idx, "GROUP");
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }
    rc = VCursorAddColumn(curs, &hdr_idx, "HDR");
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }
    rc = VCursorAddColumn(curs, &tag_idx, "TAG");
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }
    rc = VCursorAddColumn(curs, &value_idx, "VALUE");
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    rc = VCursorOpen(curs);
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }

    int64_t start = 0;
    uint64_t count = 0;
    rc = VCursorIdRange(curs, 0, &start, &count);
    if (rc) {
        fprintf(stderr, "Failed %d %d", __LINE__, rc);
        return rc;
    }
    printf("start=%ld,count=%lu\n", start, count);

    uint64_t prev_group=0;
    bool has_sqsn=false;
    while (count--) {
        uint64_t group;
        uint32_t row_len = 0;
        rc = VCursorReadDirect(curs, start, group_idx, 64, &group, 1,
                               &row_len);
        if (rc) {
            fprintf(stderr, "Failed %d %d", __LINE__, rc);
            return rc;
        }
        printf("group=%lu, row_len=%d\n", group, row_len);

        char hdr[8192];
        rc = VCursorReadDirect(curs, start, hdr_idx, 8, &hdr, 8192, &row_len);
        if (rc) {
            fprintf(stderr, "Failed %d %d", __LINE__, rc);
            return rc;
        }
        hdr[row_len] = '\0';
        printf("hdr=%s, row_len=%d\n", hdr, row_len);

        char tag[8192];
        rc = VCursorReadDirect(curs, start, tag_idx, 8, &tag, 8192, &row_len);
        if (rc) {
            fprintf(stderr, "Failed %d %d", __LINE__, rc);
            return rc;
        }
        tag[row_len] = '\0';
        printf("tag=%s, row_len=%d\n", tag, row_len);

        char value[8192];
        rc = VCursorReadDirect(curs, start, value_idx, 8, &value, 8192, &row_len);
        if (rc) {
            fprintf(stderr, "Failed %d %d", __LINE__, rc);
            return rc;
        }
        value[row_len] = '\0';
        printf("value=%s, row_len=%d\n", value, row_len);

        if (!strcmp(hdr,"SQ") && !strcmp(tag,"SN"))
        {
            if (prev_group && group==prev_group)
            {
                if (has_sqsn)
                {
                    fprintf(stderr,"Repeated SQSN\n");
                }
                has_sqsn=true;
            }
        }

        if (group!=prev_group)
        {
            prev_group=group;
            has_sqsn=false;
        }

        start++;
    }
    fprintf(stderr, "Made verifier for %s\n", dbname);

    VCursorRelease(curs);
    VTableRelease(tbl);
    VDatabaseRelease(db);
    VDBManagerRelease(mgr);
    KDirectoryRelease(srcdir);

    return 0;
}

} // namespace ncbi

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
        const char* dbname = *(++argv);
        rc = ncbi::process(dbname);
    }
    return rc;
}
}
