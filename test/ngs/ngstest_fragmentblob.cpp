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
* Unit tests for low-level NGS functions
*/

// suppress macro max from windows.h
#define NOMINMAX

#include "ngs_c_fixture.hpp"

#include <NGS_Cursor.h>
#include <SRA_Read.h>
#include <NGS_FragmentBlob.h>
#include <NGS_FragmentBlobIterator.h>

#include <vdb/table.h>
#include <vdb/database.h>

#include <stdexcept>
#include <cstring>
#include <limits>
#include <cmath>

using namespace std;
using namespace ncbi::NK;

TEST_SUITE(NgsFragmentBlobTestSuite);

const char* SRA_Accession = "SRR000001";

//////////////////////////////////////////// NGS_FragmentBlob

class FragmentBlobFixture : public NGS_C_Fixture
{
public:
    FragmentBlobFixture ()
    :   m_tbl ( 0 ),
        m_curs ( 0 ),
        m_blob ( 0 )
    {
    }

    void MakeBlob ( const char* acc, int64_t rowId )
    {
        if ( m_tbl != 0 )
            VTableRelease ( m_tbl );
        if ( VDBManagerOpenTableRead ( m_ctx -> rsrc -> vdb, & m_tbl, NULL, acc ) != 0 )
            throw logic_error ("FragmentBlobFixture::MakeCursor VDBManagerOpenTableRead failed");
        m_curs = NGS_CursorMake ( m_ctx, m_tbl, sequence_col_specs, seq_NUM_COLS );
        if ( m_curs == 0 )
            throw logic_error ("FragmentBlobFixture::MakeCursor NGS_CursorMake failed");
        NGS_String* run = NGS_StringMake ( m_ctx, acc, string_size ( acc ) );
        m_blob = NGS_FragmentBlobMake ( m_ctx, run, m_curs, rowId );
        NGS_StringRelease ( run, m_ctx );
        if ( m_blob == 0 )
            throw logic_error ("FragmentBlobFixture::MakeCursor NGS_FragmentBlobMake failed");
    }
    virtual void Release()
    {
        if (m_ctx != 0)
        {
            if ( m_blob != 0 )
            {
                NGS_FragmentBlobRelease ( m_blob, m_ctx );
            }
            if ( m_curs != 0 )
            {
                NGS_CursorRelease ( m_curs, m_ctx );
            }
            if ( m_tbl != 0 )
            {
                VTableRelease ( m_tbl );
            }
        }
        NGS_C_Fixture :: Release ();
    }

    const VTable* m_tbl;
    const NGS_Cursor* m_curs;
    struct NGS_FragmentBlob* m_blob;
};

TEST_CASE ( NGS_FragmentBlobMake_BadCursor)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );

    NGS_String* run = NGS_StringMake ( ctx, "", 0 );
    REQUIRE ( ! FAILED () );

    struct NGS_FragmentBlob * blob = NGS_FragmentBlobMake ( ctx, run, NULL, 1 );
    REQUIRE_FAILED ();
    REQUIRE_NULL ( blob );

    NGS_StringRelease ( run, ctx );
    REQUIRE ( ! FAILED () );
}

FIXTURE_TEST_CASE ( NGS_FragmentBlob_RowRange, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    int64_t first = 0;
    uint64_t count = 0;
    NGS_FragmentBlobRowRange ( m_blob, m_ctx, & first, & count );
    REQUIRE_EQ ( (int64_t)1, first );
    REQUIRE_EQ ( (uint64_t)4, count );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_BadRowId, FragmentBlobFixture )
{
    ENTRY;
    // BadRowId. The object gets created but NGS_FragmentBlobRowRange will fail
    MakeBlob ( SRA_Accession, -1 );

    int64_t first = 0;
    uint64_t count = 0;
    NGS_FragmentBlobRowRange ( m_blob, ctx, & first, & count );
    REQUIRE_FAILED ();

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_DuplicateRelease, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    // Duplicate
    NGS_FragmentBlob* anotherBlob = NGS_FragmentBlobDuplicate ( m_blob, m_ctx );
    REQUIRE_NOT_NULL ( anotherBlob );
    // Release
    NGS_FragmentBlobRelease ( anotherBlob, m_ctx );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_Data_BadArg, FragmentBlobFixture )
{
    ENTRY;

    REQUIRE_NULL ( NGS_FragmentBlobData ( NULL, m_ctx ) );
    REQUIRE_FAILED ();

    EXIT;
}
FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_Data, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    const void* data = NGS_FragmentBlobData ( m_blob, m_ctx );
    REQUIRE_NOT_NULL ( data );
    REQUIRE_EQ ( string ( "TCAGAT" ), string ( (const char*)data, 6 ) );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_Size_BadArg, FragmentBlobFixture )
{
    ENTRY;

    NGS_FragmentBlobSize ( NULL, m_ctx );
    REQUIRE_FAILED ();

    EXIT;
}
FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_Size, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    REQUIRE_EQ ( (uint64_t)1080, NGS_FragmentBlobSize ( m_blob, m_ctx ) );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_InfoByOffset_BadSelf, FragmentBlobFixture )
{
    ENTRY;

    int64_t rowId;
    uint64_t fragStart;
    uint64_t baseCount;
    int32_t bioNumber;
    NGS_FragmentBlobInfoByOffset ( NULL, ctx, 0, & rowId, & fragStart, & baseCount, & bioNumber );
    REQUIRE_FAILED ();

    EXIT;
}
// TODO: NULL for optional parameters

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_InfoByOffset_Biological, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    int64_t rowId;
    uint64_t fragStart;
    uint64_t baseCount;
    int32_t bioNumber;
    // offset 300 is in row 2 which starts at 284 and consists of 4 fragments:
    // technical, start 284, len 4
    // biological #0, start 288, len 115 <== expect to see this for offset 300
    // technical, start 403, len 44
    // biological #1, start 447, len 99
    NGS_FragmentBlobInfoByOffset ( m_blob, ctx, 300, & rowId, & fragStart, & baseCount, & bioNumber );
    REQUIRE_EQ ( (int64_t)2, rowId );
    REQUIRE_EQ ( (uint64_t)288, fragStart );
    REQUIRE_EQ ( (uint64_t)115, baseCount );
    REQUIRE_EQ ( (int32_t)0, bioNumber );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobMake_InfoByOffset_Technical, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    int64_t rowId;
    uint64_t fragStart;
    uint64_t baseCount;
    int32_t bioNumber;
    // offset 300 is in row 2 which starts at 284 and consists of 4 fragments:
    // technical, start 284, len 4
    // biological #0, start 288, len 115
    // technical, start 403, len 44  <== expect to see this for offset 410
    // biological #1, start 447, len 99
    NGS_FragmentBlobInfoByOffset ( m_blob, ctx, 410, & rowId, & fragStart, & baseCount, & bioNumber );
    REQUIRE_EQ ( (int64_t)2, rowId );
    REQUIRE_EQ ( (uint64_t)403, fragStart );
    REQUIRE_EQ ( (uint64_t)44, baseCount );
    REQUIRE_EQ ( (int32_t)-1, bioNumber );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlob_MakeFragmentId, FragmentBlobFixture )
{
    ENTRY;
    MakeBlob ( SRA_Accession, 1 );

    NGS_String* id = NGS_FragmentBlobMakeFragmentId ( m_blob, ctx, 2, 1 );
    REQUIRE_EQ ( string ( SRA_Accession ) + ".FR1.2", string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );
    NGS_StringRelease ( id, ctx );

    EXIT;
}

//////////////////////////////////////////// NGS_FragmentBlobIterator

class BlobIteratorFixture : public NGS_C_Fixture
{
public:
    BlobIteratorFixture ()
    :   m_tbl ( 0 ),
        m_blobIt ( 0 )
    {
    }

    void MakeSRA( const char* acc )
    {
        if ( m_tbl != 0 )
            VTableRelease ( m_tbl );
        if ( VDBManagerOpenTableRead ( m_ctx -> rsrc -> vdb, & m_tbl, NULL, acc ) != 0 )
            throw logic_error ("BlobIteratorFixture::MakeSRA VDBManagerOpenTableRead failed");
    }
    void MakeIterator( const char* acc )
    {
        MakeSRA ( acc );
        NGS_String* run = NGS_StringMake ( m_ctx, acc, string_size ( acc ) );
        m_blobIt = NGS_FragmentBlobIteratorMake ( m_ctx, run, m_tbl );
        NGS_StringRelease ( run, m_ctx );
    }
    virtual void Release()
    {
        if (m_ctx != 0)
        {
            if ( m_blobIt != 0 )
            {
                NGS_FragmentBlobIteratorRelease ( m_blobIt, m_ctx );
            }
            if ( m_tbl != 0 )
            {
                VTableRelease ( m_tbl );
            }
        }
        NGS_C_Fixture :: Release ();
    }

    const VTable* m_tbl;
    struct NGS_FragmentBlobIterator* m_blobIt;
};

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_BadMake, BlobIteratorFixture )
{
    ENTRY;

    NGS_String* run = NGS_StringMake ( m_ctx, "", 0 );
    REQUIRE ( ! FAILED () );

    struct NGS_FragmentBlobIterator* blobIt = NGS_FragmentBlobIteratorMake ( ctx, run, NULL );
    REQUIRE_FAILED ();
    REQUIRE_NULL ( blobIt );
    NGS_StringRelease ( run, ctx );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_CreateRelease, BlobIteratorFixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession );

    NGS_String* run = NGS_StringMake ( m_ctx, SRA_Accession, string_size ( SRA_Accession ) );
    REQUIRE ( ! FAILED () );
    struct NGS_FragmentBlobIterator* blobIt = NGS_FragmentBlobIteratorMake ( m_ctx, run, m_tbl );
    REQUIRE ( ! FAILED () );
    NGS_StringRelease ( run, m_ctx );
    REQUIRE_NOT_NULL ( blobIt );
    REQUIRE ( ! FAILED () );
    // Release
    NGS_FragmentBlobIteratorRelease ( blobIt, m_ctx );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_DuplicateRelease, BlobIteratorFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession );

    // Duplicate
    struct NGS_FragmentBlobIterator* anotherBlobIt = NGS_FragmentBlobIteratorDuplicate ( m_blobIt, m_ctx );
    REQUIRE_NOT_NULL ( anotherBlobIt );
    // Release
    NGS_FragmentBlobIteratorRelease ( anotherBlobIt, m_ctx );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_Next, BlobIteratorFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession );

    struct NGS_FragmentBlob* blob = NGS_FragmentBlobIteratorNext ( m_blobIt, m_ctx );
    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( blob );
    NGS_FragmentBlobRelease ( blob, ctx );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_HasMore, BlobIteratorFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession );

    REQUIRE ( NGS_FragmentBlobIteratorHasMore ( m_blobIt, m_ctx ) );
    REQUIRE ( ! FAILED () );

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_FragmentBlobIterator_FullScan, BlobIteratorFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession );

    uint32_t count = 0;
    while ( NGS_FragmentBlobIteratorHasMore ( m_blobIt, m_ctx ) )
    {
        struct NGS_FragmentBlob* blob = NGS_FragmentBlobIteratorNext ( m_blobIt, m_ctx );
        REQUIRE_NOT_NULL ( blob );
        NGS_FragmentBlobRelease ( blob, ctx );
        ++count;
    }
    REQUIRE_EQ ( (uint32_t)243, count);
    REQUIRE_NULL ( NGS_FragmentBlobIteratorNext ( m_blobIt, m_ctx ) );

    EXIT;
}


//////////////////////////////////////////// Main

extern "C"
{

#include <kapp/args.h>

ver_t CC KAppVersion ( void )
{
    return 0x1000000;
}
rc_t CC UsageSummary (const char * progname)
{
    return 0;
}

rc_t CC Usage ( const Args * args )
{
    return 0;
}

const char UsageDefaultName[] = "test-ngs";

rc_t CC KMain ( int argc, char *argv [] )
{
    rc_t m_coll=NgsFragmentBlobTestSuite(argc, argv);
    return m_coll;
}

}

