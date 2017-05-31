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

#include <../../libs/align/samextract-pool.h>
#include <../../libs/align/samextract.h>
#include <align/samextract-lib.h>
#include <ctype.h>
#include <kapp/args.h>
#include <kapp/main.h>
#include <kfg/config.h>
#include <kfs/directory.h>
#include <kfs/file.h>
#include <klib/defs.h>
#include <klib/rc.h>
#include <klib/text.h>
#include <klib/vector.h>
#include <ktst/unit_test.hpp>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strtol.h>
#include <time.h>
#include <unistd.h>

#include <cstdlib>
#include <stdexcept>
#include <sysalloc.h>

using namespace std;

static const size_t NUM_RAND = 100000;
//#define TEST_ALL_THE_INTEGERS // Takes about 2.5 hours

static bool tst_fast_u32toa(u32 val)
{
    char slow[100];
    char fast[100];
    sprintf(slow, "%u", val);
    fast_u32toa(fast, val);
    if (strcmp(fast, slow)) {
        fprintf(stderr, "mismatch %u '%s' '%s'\n", val, slow, fast);
        return false;
    }
    return true;
}

static bool tst_fast_i32toa(u32 val)
{
    char slow[100];
    char fast[100];
    sprintf(slow, "%d", val);
    fast_i32toa(fast, val);
    if (strcmp(fast, slow)) {
        fprintf(stderr, "mismatch %d '%s' '%s'\n", val, slow, fast);
        return false;
    }
    return true;
}

TEST_SUITE(SAMExtractTestSuite)

TEST_CASE(Fast_Sequence)
{
    unsigned char test[]
        = {0x00, 0xFF, 0x01, 0x10, 0x0F, 0xF0, 0x00, 0x01, 0x01, 0x01};
    char dest[30];

    decode_seq(test, 0, dest);
    REQUIRE_EQUAL(strcmp(dest, ""), 0);
    decode_seq(test, 1, dest);
    REQUIRE_EQUAL(strcmp(dest, "="), 0);
    decode_seq(test, 2, dest);
    REQUIRE_EQUAL(strcmp(dest, "=="), 0);
    decode_seq(test, 3, dest);
    REQUIRE_EQUAL(strcmp(dest, "==N"), 0);
    decode_seq(test, 4, dest);
    REQUIRE_EQUAL(strcmp(dest, "==NN"), 0);
    decode_seq(test, 5, dest);
    REQUIRE_EQUAL(strcmp(dest, "==NN="), 0);
    decode_seq(test, 20, dest);
    REQUIRE_EQUAL(strcmp(dest, "==NN=AA==NN====A=A=A"), 0);
}

TEST_CASE(Fast_u32toa)
{
    static const u32 tsts[]
        = {0,          1,          2,           9,           10,
           11,         99,         100,         101,         110,
           199,        200,        201,         999,         1000,
           1001,       1010,       1019,        1900,        1988,
           1990,       1992,       2003,        2005,        9999,
           10000,      99999,      20030613,    20050722,    100000,
           1000000,    10000000,   10000000,    100000000,   1000000000,
           2000000000, 2147483647, 2147483648l, 4294967294l, 4294967295l};
    for (int i = 0; i != sizeof(tsts) / sizeof(tsts[0]); ++i)
        REQUIRE_EQUAL(tst_fast_u32toa(tsts[i]), true);

    for (int i = 0; i != NUM_RAND; ++i)
        REQUIRE_EQUAL(tst_fast_u32toa(lrand48()), true);

#ifdef TEST_ALL_THE_INTEGERS
    fprintf(stderr, "all u32toa\n");
    for (u32 u = 0; u != UINT32_MAX; ++u)
        REQUIRE_EQUAL(tst_fast_u32toa(u), true);
#endif
}

static bool tst_strtoi64(const char* str)
{
    i64 slow = strtoi64(str, NULL, 10);
    i64 fast = fast_strtoi64(str);
    if (fast != slow) {
        fprintf(stderr, "mismatch '%s' %ld %ld\n", str, slow, fast);
        return false;
    }
    return true;
}

TEST_CASE(Fast_i32toa)
{
    static const i32 tsts[]
        = {0,         1,          2,          9,          10,        11,
           99,        199,        200,        999,        1000,      9999,
           10000,     99999,      100000,     1000000,    10000000,  10000000,
           100000000, 1000000000, 2000000000, 2147483646, 2147483647};

    for (int i = 0; i != sizeof(tsts) / sizeof(tsts[0]); ++i) {
        REQUIRE_EQUAL(tst_fast_i32toa(tsts[i]), true);
        REQUIRE_EQUAL(tst_fast_i32toa(-tsts[i]), true);
    }

    for (int i = 0; i != NUM_RAND; ++i)
        REQUIRE_EQUAL(tst_fast_i32toa(mrand48()), true);

#ifdef TEST_ALL_THE_INTEGERS
    fprintf(stderr, "all i32toa\n");
    for (int i = INT32_MIN; i != INT32_MAX; ++i)
        REQUIRE_EQUAL(tst_fast_i32toa(i), true);
#endif
}

TEST_CASE(Fast_strtoi64)
{
    static const char* tsts[]
        = {"0",          "1",          "9",          "10",
           "99",         "100",        "101",        "199",
           "999",        "999999",     "1000000",    "100000000",
           "2147483645", "2147483646", "8589934592", "9223372036854775807"};
    char str[32];

    for (int i = 0; i != sizeof(tsts) / sizeof(tsts[0]); ++i) {
        REQUIRE_EQUAL(tst_strtoi64(tsts[i]), true);
        sprintf(str, "-%s", tsts[i]);
        REQUIRE_EQUAL(tst_strtoi64(str), true);
    }

    for (int i = 0; i != NUM_RAND; ++i) {
        sprintf(str, "%ld", mrand48() * lrand48() + mrand48());
        REQUIRE_EQUAL(tst_strtoi64(str), true);
    }
#ifdef TEST_ALL_THE_INTEGERS
    fprintf(stderr, "all strtoi64\n");
    for (int i = INT32_MIN; i != INT32_MAX; ++i) {
        sprintf(str, "%d", i);
        REQUIRE_EQUAL(tst_strtoi64(str), true);
    }
#endif
}

TEST_CASE(Decode_Cigar)
{
    pool_init();
    u32   incigar1[] = {0x10, 0x21, 0x38};
    char* outcigar; // pool allocated
    outcigar = decode_cigar(incigar1, sizeof(incigar1) / sizeof(incigar1[0]));
    REQUIRE_EQUAL(strcmp("1M2I3X", outcigar), 0);
    u32 incigar2[] = {0xfffffff0};
    outcigar = decode_cigar(incigar2, sizeof(incigar2) / sizeof(incigar2[0]));
    REQUIRE_EQUAL(strcmp("268435455M", outcigar), 0);
    pool_release();
}

TEST_CASE(In_Range)
{
    REQUIRE_EQUAL(inrange("0", 0, 0), true);
    REQUIRE_EQUAL(inrange("0", 1, 0), false);
    REQUIRE_EQUAL(inrange("0", 1, 2), false);
    REQUIRE_EQUAL(inrange("0", -1, -1), false);
}

TEST_CASE(Is_MD5)
{
    REQUIRE_EQUAL(ismd5(""), false);
    REQUIRE_EQUAL(ismd5("A"), false);
    REQUIRE_EQUAL(ismd5("-"), false);
    REQUIRE_EQUAL(ismd5("7-ce4e3aa8b07d1f448ace526a3977fd"), false);
    REQUIRE_EQUAL(ismd5("7fce4e3aa8b07d1f448ace526a3977fd"), true);
    REQUIRE_EQUAL(ismd5("7fce4e3aa8b07d1f448ace526a3977f*"), true);
}

TEST_CASE(Is_floworder)
{
    REQUIRE_EQUAL(isfloworder("*"), true);
    REQUIRE_EQUAL(isfloworder("ACMGRSVTWYHKDBN"), true);
    REQUIRE_EQUAL(isfloworder("0"), false);
    REQUIRE_EQUAL(isfloworder("a"), false);
}
// filter
// mempool

extern "C" {
ver_t CC KAppVersion(void) { return 0x1000000; }
rc_t CC UsageSummary(const char* progname) { return 0; }

rc_t CC Usage(const Args* args) { return 0; }

const char UsageDefaultName[] = "test-samextract";

rc_t CC KMain(int argc, char* argv[])
{
    srand48(time(NULL));
    rc_t rc = SAMExtractTestSuite(argc, argv);
    return rc;
}
}
