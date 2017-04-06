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

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <klib/vector.h>
#include "samextract.h"
#include <stdint.h>

const size_t BLOCK_SZ=2 * 1024 * 1024;

Vector allocs;

void * cur_block=NULL;
size_t cur_block_remain;

static void morecore(size_t sz)
{
    cur_block=malloc(sz);
    if (cur_block==NULL)
    {
        ERR("out of memory");
    }
    cur_block_remain=sz;

    VectorAppend(&allocs,NULL,cur_block);
    DBG("morecore %p", cur_block);
}

void pool_init(void)
{
    DBG("pool_init");
    if (cur_block)
        ERR("Need to release pool");

    VectorInit(&allocs,0,0);
    morecore(BLOCK_SZ);
}

void pool_release(void)
{
    DBG("pool_release");
    DBG("pool used %d", BLOCK_SZ-cur_block_remain);
    DBG("pools:%d", VectorLength(&allocs));
    for (u32 i=0; i!=VectorLength(&allocs); ++i)
    {
        DBG("freeing %p", VectorGet(&allocs,i));
        free(VectorGet(&allocs,i));
    }

    VectorWhack(&allocs,NULL,NULL);
    cur_block=NULL;
    pool_init();
}

void pool_free(void * buf)
{
    /* Could conceivably reclaim last allocation */
    return;
}

void * pool_alloc(size_t sz)
{
    if (sz % 8 != 0) sz += 8-(sz % 8); /* Round up for alignment */
    if (sz > cur_block_remain)
        morecore(MAX(sz,BLOCK_SZ));

    void * buf=cur_block;
    cur_block += sz;
    cur_block_remain -= sz;

    return buf;
}

char * pool_strdup(const char * str)
{
    size_t len=strlen(str)+1;
    void * buf=pool_alloc(len);
    memmove(buf,str,len);
    return (char *)buf;
}

char * pool_memdup(const char * str, size_t len)
{
    void * buf=pool_alloc(len);
    memmove(buf,str,len);
    return (char *)buf;
}

