// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//               National Center for Biotechnology Information
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act.  It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted.  This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government have not placed any restriction on its use or reproduction.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
//  Please cite the author in any work or product based on this material.
//
// ===========================================================================

#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/file.h>
#include <kfs/directory.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <klib/text.h>
#include <align/samextract-lib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <ktst/unit_test.hpp> // TEST_CASE
#include <kfg/config.h>

#include <sysalloc.h>
#include <cstdlib>
#include <stdexcept>

using namespace std;

TEST_SUITE(SAMExtractTestSuite)

TEST_CASE(SAMExtract)
{
}

/*
 * fast sequence extract tests:
    // Do some self-checks during initialization
    if (memcmp(&seqbytemap[0], "==", 2) || memcmp(&seqbytemap[1], "=A", 2)
        || memcmp(&seqbytemap[16], "A=", 2)
        || memcmp(&seqbytemap[255], "NN", 2))
    {
        ERR("Self-check failed: seqbytemap %s", seqbytemap);
    }

    //           ==    NN    =A    A=    =N    N=    ==    =A    =A    =A
    u8 test[] = {0x00, 0xFF, 0x01, 0x10, 0x0F, 0xF0, 0x00, 0x01, 0x01, 0x01};
    char dest[30];

    decode_seq(test, 0, dest);
    if (strcmp(dest, "")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 1, dest);
    if (strcmp(dest, "=")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 2, dest);
    if (strcmp(dest, "==")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 3, dest);
    if (strcmp(dest, "==N")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 4, dest);
    if (strcmp(dest, "==NN")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 5, dest);
    if (strcmp(dest, "==NN=")) ERR("Self-check failed: %s", dest);
    decode_seq(test, 20, dest);
    if (strcmp(dest, "==NN=AA==NN====A=A=A"))
        ERR("Self-check failed: %s", dest);

    DBG("Self-checks OK");

*/
// mempool tests
// sequence decoding
// cigar decoding
extern "C"
{
ver_t CC KAppVersion(void) { return 0x1000000; }
rc_t CC UsageSummary(const char* progname) { return 0; }

rc_t CC Usage(const Args* args) { return 0; }

const char UsageDefaultName[] = "test-samextract";

rc_t CC KMain(int argc, char* argv[])
{
    rc_t rc = SAMExtractTestSuite(argc, argv);
    return rc;
}
}
