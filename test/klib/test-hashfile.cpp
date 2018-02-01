/*===========================================================================
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
/**
* Unit tests for hash files
*/

#include <ktst/unit_test.hpp>

#include <kfs/defs.h>
#include <kfs/directory.h>
#include <kfs/file.h>
#include <klib/data-buffer.h>
#include <klib/hashfile.h>
#include <klib/hashtable.h>
#include <klib/log.h>
#include <klib/misc.h> /* is_user_admin() */
#include <klib/num-gen.h>
#include <klib/printf.h>
#include <klib/sort.h>
#include <klib/text.h>
#include <klib/vector.h>

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/time.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

static const size_t RAND_SIZE = 1048576;
static char* RANDS = NULL;

/* #define BENCHMARK */

TEST_SUITE(KHashFileTestSuite);

TEST_CASE(Klib_KHashFileSet)
{
    rc_t rc;
    const char* str1 = "Tu estas probando este hoy, no manana";
    const char* str2 = "Tu estas probando este hoy, no mananX";
    size_t size = strlen(str1);

    KHashFile* hset = NULL;

    rc = KHashFileMake(NULL, NULL);
    REQUIRE_RC_FAIL(rc);
    REQUIRE_EQ((void*)hset, (void*)NULL);

    rc = KHashFileMake(&hset, NULL);
    REQUIRE_RC(rc);
    REQUIRE_NE((void*)hset, (void*)NULL);

    size_t sz = KHashFileCount(hset);
    REQUIRE_EQ(sz, (size_t)0);

    uint64_t hash = KHash(str1, size);
    rc = KHashFileAdd(hset, str1, strlen(str1), hash, NULL, 0);
    REQUIRE_RC(rc);

    sz = KHashFileCount(hset);
    REQUIRE_EQ(sz, (size_t)1);

    rc = KHashFileAdd(hset, str1, strlen(str1), hash, NULL, 0);
    REQUIRE_RC(rc);

    sz = KHashFileCount(hset);
    REQUIRE_EQ(sz, (size_t)1);

    bool found;
    found = KHashFileFind(hset, str1, strlen(str1), hash, NULL, NULL);
    REQUIRE_EQ(found, true);

    found = KHashFileFind(hset, str2, strlen(str2), hash, NULL, NULL);
    REQUIRE_EQ(found, false);

    KHashFileDispose(hset);
}

TEST_CASE(Klib_hashfileMap)
{
    const char* str1 = "Tu estas probando este hoy, no manana";
    const char* str2 = "Tu estas probando este hoy, no mananX";

    KDirectory* dir = NULL;
    rc_t rc;
    rc = KDirectoryNativeDir(&dir);
    REQUIRE_RC(rc);

    const char* fname = tmpnam(NULL);
    KFile* backing = NULL;
    rc = KDirectoryCreateFile(dir, &backing, true, 0600, kcmInit, fname);
    // TODO | kcmParents?
    REQUIRE_RC(rc);
//    rc = KFileSetSize(backing, 4096);
//    REQUIRE_RC(rc);

    KHashFile* hmap;
    rc = KHashFileMake(&hmap, backing);
    REQUIRE_RC(rc);

    size_t sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)0);

    uint64_t hash = 1;
    uint64_t val1 = 123;
    rc = KHashFileAdd(hmap, str1, strlen(str1), hash, &val1, sizeof(val1));
    REQUIRE_RC(rc);

    sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)1);

    rc = KHashFileAdd(hmap, str1, strlen(str1), hash, &val1, sizeof(val1));
    REQUIRE_RC(rc);

    sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)1);

    bool found;
    uint64_t val;
    uint64_t len = 0;
    found = KHashFileFind(hmap, str1, strlen(str1), hash, NULL, 0);
    REQUIRE_EQ(found, true);
    found = KHashFileFind(hmap, str1, strlen(str1), hash, &val, &len);
    REQUIRE_EQ(found, true);
    REQUIRE_EQ(len, sizeof(val));
    len = 0;
    REQUIRE_EQ(val, (uint64_t)123);

    uint64_t val2 = 124;
    rc = KHashFileAdd(hmap, str1, strlen(str1), hash, &val2, sizeof(val2));
    REQUIRE_RC(rc);

    sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)1);

    found = KHashFileFind(hmap, str1, strlen(str1), hash, &val, &len);
    REQUIRE_EQ(found, true);
    REQUIRE_EQ(len, sizeof(val));
    len = 0;
    REQUIRE_EQ(val, (uint64_t)124);

    found = KHashFileFind(hmap, str2, strlen(str2), hash, NULL, 0);
    REQUIRE_EQ(found, false);

    uint64_t val3 = 125;
    rc = KHashFileAdd(hmap, str2, strlen(str2), hash, &val3, sizeof(val3));
    REQUIRE_RC(rc);

    found = KHashFileFind(hmap, str2, strlen(str2), hash, &val, &len);
    REQUIRE_EQ(found, true);
    REQUIRE_EQ(len, sizeof(val));
    len = 0;
    REQUIRE_EQ(val, (uint64_t)125);

    sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)2);

    KHashFileDispose(hmap);
    rc = KDirectoryRemove(dir, true, "%s", fname);
    REQUIRE_RC(rc);
    KFileRelease(backing);
    KDirectoryRelease(dir);
}

TEST_CASE(Klib_hashfileMapDeletes)
{
    rc_t rc;

    KDirectory* dir = NULL;
    rc = KDirectoryNativeDir(&dir);
    REQUIRE_RC(rc);

    const char* fname = tmpnam(NULL);
    KFile* backing = NULL;
    rc = KDirectoryCreateFile(dir, &backing, true, 0600, kcmInit, fname);
    REQUIRE_RC(rc);

    KHashFile* hmap = NULL;
    rc = KHashFileMake(&hmap, backing);
    REQUIRE_RC(rc);

    size_t sz = KHashFileCount(hmap);
    REQUIRE_EQ(sz, (size_t)0);

    std::vector<std::string> strs, vals;
    const size_t loops = 100000;
    for (size_t i = 0; i != loops; ++i) {
        strs.push_back(
            string(1 + (random() % 500), char(32 + random() % 90)));
        vals.push_back(
            string(1 + (random() % 500), char(32 + random() % 90)));
    }

    std::unordered_map<std::string, std::string> map;
    for (size_t i = 0; i != strs.size(); ++i) {
        const auto key = strs[i];
        const auto val = vals[i];
        auto pair = std::make_pair(key, val);
        map.erase(key);
        map.insert(pair);

        uint64_t hash = KHash(key.data(), key.size());
        KHashFileDelete(hmap, key.data(), key.size(), hash);
        rc = KHashFileAdd(hmap, key.data(), key.size(), hash, val.data(),
                          val.size());
        REQUIRE_RC(rc);

        if (random() % 2) {
            map.erase(key);
            bool found = KHashFileDelete(hmap, key.data(), key.size(), hash);
            REQUIRE_EQ(found, true);
        }

        sz = KHashFileCount(hmap);
        REQUIRE_EQ(sz, (size_t)map.size());
    }

    for (auto it : map) {
        const auto key = it.first;
        const auto value = it.second;

        uint64_t hash = KHash(key.data(), key.size());
        char val[500];
        size_t len = 0;
        bool found
            = KHashFileFind(hmap, key.data(), key.size(), hash, &val, &len);
        REQUIRE_EQ(found, true);
        REQUIRE_EQ(value.size(), len);
        REQUIRE_EQ((int)0, memcmp(value.data(), val, len));
    }

    KHashFileDispose(hmap);
    rc = KDirectoryRemove(dir, true, "%s", fname);
    REQUIRE_RC(rc);
    KFileRelease(backing);
    KDirectoryRelease(dir);
}

TEST_CASE(Klib_hashfilethreads)
{

}

#ifdef BENCHMARK
#endif // BENCHMARK

extern "C" {

#include <kapp/args.h>
#include <kfg/config.h>

ver_t CC KAppVersion(void) { return 0x1000000; }
rc_t CC UsageSummary(const char* progname) { return 0; }

rc_t CC Usage(const Args* args) { return 0; }

const char UsageDefaultName[] = "test-hashfile";

rc_t CC KMain(int argc, char* argv[])
{
    srandom(time(NULL));

    RANDS = (char*)malloc(RAND_SIZE);

    for (size_t i = 0; i != RAND_SIZE; ++i) RANDS[i] = random();

    KConfigDisableUserSettings();
    rc_t rc = KHashFileTestSuite(argc, argv);

    free(RANDS);
    return rc;
}
}
