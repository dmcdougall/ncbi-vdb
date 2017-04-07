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

#define __STDC_LIMIT_MACROS
#include <kapp/args.h>
#include <kapp/main.h>
#include <kfs/file.h>
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
#include <strtol.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include "samextract.h"
#include "samextract-tokens.h"
#include <align/samextract-lib.h>
#include <stdint.h>

// TODO: Put these into struct
static const ssize_t READBUF_SZ=65536;
static char readbuf[READBUF_SZ+1];
static size_t readbuf_sz=0;
static size_t readbuf_pos=0;

static char curline[READBUF_SZ+1];
static size_t curline_len=0;
static uint64_t file_pos=0;

class chunkview
{
  public:
      chunkview(KQueue * queue) : que(queue), c(NULL), cur(NULL) {};

      ~chunkview()
      {
          que=NULL;
          c=NULL;
          cur=NULL;
      }

      // TODO: Handle _MSC_VER
      // TODO: Move to a USE_ include that defines features:
      // __HAS_RVALUE_REFERENCES, ...
      // __HAS_METHOD_DELETE, ...
      // __has_cpp_attribute
      // __CPPVER=98,11,14,..
#ifndef __cplusplus
      // C++98
      chunkview(const chunkview &);
#elif __cplusplus <= 199711L
      // C++98
      chunkview(const chunkview &);
#elif __cplusplus >= 201103L
      // C++11, C++14
      chunkview(const chunkview &) = delete; // No copy ctor
      chunkview & operator=(const chunkview &) = delete; // No assignment
      chunkview(const chunkview &&) = delete; // No move ctor
      chunkview & operator=(const chunkview &&) = delete; // No move assignment
#endif

  private:
      bool getnextchunk()
      {
          void * where=NULL;
          struct timeout_t tm;
          DBG("\t\tBlocker thread checking blocker queue");
          if (c && c->state!=uncompressed)
          {
              ERR("\t\tblocker bad state");
              return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
          }
          if (c) c->state=empty;
          c=NULL;

          while (1)
          {
              TimeoutInit(&tm, 5000); // 5 seconds
              rc_t rc=KQueuePop(que, &where, &tm);
              if (rc==0)
              {
                  c=(chunk *)where;
                  if (c->state==empty)
                  {
                      ERR("\t\tBlocker bad state");
                      return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
                  }
                  while (c->state!=uncompressed)
                  {
                      DBG("\t\tBlocker busy");
                      usleep(1);
                  }
                  DBG("\t\tBlocker thread chunk %p size %u", c->in, c->insize);
                  break;
              } else if ((int)GetRCObject(rc)== rcTimeout)
              {
                  INFO("\t\tBlocker thread queue empty");
                  if (KQueueSealed(que))
                  {
                      INFO("\t\tQueue sealed, Blocker thread complete");
                      return false;
                  }
                  //            usleep(1);
              } else if ((int)GetRCObject(rc)== rcData)
              {
                  INFO("\t\tBlocker thread queue data (empty?)");
                  return false;
              } else
              {
                  ERR("blocker rc=%d",rc);
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
              }
          }
          return true;
      }
  public:
      bool getbytes(void * dest, u8 len)
      {
          char * where=(char *)dest;
          while (len)
          {
              if (c==NULL || c->outsize==0)
              {
                  if (!getnextchunk())
                      return false;
                  cur=c->out;
              }

              size_t howmany=c->outsize;
              if (len < c->outsize) howmany=len;

              memmove(where,cur,howmany);
              len-=howmany;
              c->outsize-=howmany;
              where+=howmany;
              cur+=howmany;
          }
          return true;
      }

  private:
      KQueue * que;
      chunk * c;
      Bytef * cur;
};


void logmsg (const char * fname, int line, const char * func, const char * severity, const char * fmt, ...)
{
    char * buf;
    size_t bufsize=0;
    FILE * buffd;

    const char * basename=strrchr(fname,'/');
    if (!basename) basename=strrchr(fname,'\\');
    if (basename) ++basename;
    if (!basename) basename=fname;
    va_list args;
    va_start(args, fmt);

    buffd=open_memstream(&buf,&bufsize);
    if (buffd==NULL)
    {
        fprintf(stderr,"can't open memstream\n");
        return;
    }
    fprintf(buffd, "%s:", severity);
    vfprintf(buffd, fmt, args);
    va_end(args);
    fprintf(buffd, "\t[%s:%s():%d]\n", basename, func, line);
    fclose(buffd);
    size_t r=fwrite(buf, bufsize, 1, stderr);
    if (r!=1) fprintf(stderr,"previous %zd log message truncated\n", bufsize);
    free(buf);
    buf=NULL;
    fflush(stderr);
    if (!strcmp(severity,"Error")) abort();
}

rc_t SAM_parseline(Extractor * state)
{
    state->rc=0;
    DBG("Parsing line (%d bytes): '%s'", strlen(curline), curline);
    SAMparse(state);
    return state->rc;
}

int moredata(char * buf, int * numbytes, size_t maxbytes)
{
    if (!curline_len)
        DBG("nomoredata");
    else
        DBG("  moredata %p %d\ncurline:'%s'", buf, maxbytes, curline);
    memmove(buf,curline,curline_len);
    *numbytes=curline_len;
    curline_len=0;
    return 0;
}

static inline rc_t readfile(Extractor * state)
{
    if (readbuf_pos == readbuf_sz)
    {
        readbuf_sz=READBUF_SZ;
        DBG("reading in at %d",file_pos);
        rc_t rc=KFileReadAll(state->infile, file_pos, readbuf, readbuf_sz, &readbuf_sz);
        file_pos+=readbuf_sz;
        if (rc) { ERR("readfile error"); return rc;}
        DBG("Read in %d",readbuf_sz);
        readbuf_pos=0;
        if (!readbuf_sz)
        {
            DBG("Buffer complete. EOF");
        }
    }
    return 0;
}

void SAMerror(Extractor * state, const char * s)
{
    ERR(" Parsing error: %s\nLine was:'%s'",s, curline);
    rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
    state->rc=rc;
    abort();
}

/* low<=str<=high */
bool inrange(const char * str, i64 low, i64 high)
{
    i64 i=strtoi64(str, NULL, 10);
    if (errno) return false;
    if (i<low || i>high) return false;
    return true;
}

bool ismd5(const char * str)
{
    size_t i;
    size_t len=strlen(str);

    if (len!=32) return false;

    for (i=0; i!=len; ++i)
    {
        if (!isalnum(str[i]) && str[i]!='*') return false;
    }

    return true;
}

/* Avoiding handling this as flex token because so many other ReadGroup values
 * could end up looking like flow orders.
 */
bool isfloworder(const char * str)
{
    size_t i;
    size_t len=strlen(str);

    if (len==1 && str[0]=='*') return true;
    for (i=0; i!=len; ++i)
    {
        switch(str[i])
        {
          case 'A':
          case 'C':
          case 'M':
          case 'G':
          case 'R':
          case 'S':
          case 'V':
          case 'T':
          case 'W':
          case 'Y':
          case 'H':
          case 'K':
          case 'D':
          case 'B':
          case 'N':
              continue;
          default:
              return false;
        }
    }
    return true;
}

rc_t process_header(Extractor * state, const char * type, const char * tag, const char * value)
{
    DBG("processing type:%s tag:%s value:%s", type, tag, value);
    if (strcmp(type,"HD") &&
        strcmp(type,"SQ") &&
        strcmp(type,"RG") &&
        strcmp(type,"PG"))
    {
        ERR("record '%s' must be HD, SQ, RG or PG", type);
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
    }

    if (strlen(tag)!=2)
    {
        ERR("tag '%s' must be 2 characters", tag);
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
    }

    if (islower(tag[0] &&
                islower(tag[1])))
    {
        DBG("optional tag");
    }

    TagValue * tv=(TagValue *)pool_alloc(sizeof(TagValue));
    if (tv==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    tv->tag=pool_strdup(tag);
    if (tv->tag==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    tv->value=pool_strdup(value);
    if (tv->value==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    VectorAppend(&state->tagvalues,NULL,tv);

    return 0;
}

rc_t mark_headers(Extractor * state, const char * type)
{
    DBG("mark_headers");
    Header * hdr=(Header *)pool_alloc(sizeof(Header));
    if (hdr==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    hdr->headercode=type;
    VectorCopy(&state->tagvalues,&hdr->tagvalues);
    VectorAppend(&state->headers,NULL,hdr);
    VectorWhack(&state->tagvalues,NULL,NULL);
    return 0;
}

rc_t process_alignment(Extractor * state, const char * qname,const char * flag,const
                       char * rname,const char * pos,const char * mapq,const char * cigar,const char *
                       rnext,const char * pnext,const char * tlen,const char * seq,const char * qual)
{
    Alignment * align=(Alignment *)pool_alloc(sizeof(Alignment));
    if (align==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
    }

    DBG("process_alignment %s %s %s", qname, rnext, qual); // TODO silence warning for now

    if (!inrange(flag,0,UINT16_MAX))
        ERR("Flag not in range %s",flag);

    if (!inrange(pos,0,INT32_MAX))
        ERR("POS not in range %s",pos);

    if (!inrange(mapq,0,UINT8_MAX))
        ERR("MAPQ not in range %s",mapq);

    if (!inrange(pnext,0,INT32_MAX))
        ERR("PNEXT not in range %s", pnext);

    if (!inrange(tlen,INT32_MIN,INT32_MAX))
        ERR("TLEN not in range %s", tlen);

    align->read=pool_strdup(seq);
    align->cigar=pool_strdup(cigar);
    align->rname=pool_strdup(rname);
    align->pos=strtou32(pos,NULL,10);
    align->flags=strtou32(flag,NULL,10);
    VectorAppend(&state->alignments,NULL,align);

    return 0;
}

// Reads next line into curline, returns false if file complete.
static bool readline(Extractor * s)
{
    DBG("readline");
    if (readfile(s)) return false;
    char * line=curline;
    line[0]='\0';
    curline_len=0;
    // Is there a newline in current buffer?
    char * nl=(char *)memchr( (readbuf+readbuf_pos), '\n', (readbuf_sz-readbuf_pos) );
    if (nl)
    {
        nl+=1;
        size_t len=nl - (readbuf+readbuf_pos);
        memmove(line, readbuf+readbuf_pos, len);
        curline_len+=len;
        readbuf_pos+=len;
        line[curline_len+1]='\0';
        return true;
    }

    // Nope, append and get more
    size_t len=(readbuf_sz - readbuf_pos);
    DBG("readline more %d/%d", readbuf_pos, len);
    memmove(line, readbuf+readbuf_pos, len);
    DBG("moreline was %d '%s'", strlen(line), line);
    line+=len;
    curline_len+=len;

    readbuf_pos=readbuf_sz;
    if (readfile(s)) return false;

    // Better be a newline now
    nl=(char *)memchr(readbuf, '\n', readbuf_sz);
    if (!nl)
    {
        return false;
    }
    DBG("found newline at %d", nl-readbuf);
    nl+=1;
    len=(nl - readbuf);
    memmove(line, readbuf, len);
    curline_len+=len;
    readbuf_pos+=len;
    line[curline_len+1]='\0';
    DBG("moreline  is %d %d '%s'", curline_len, strlen(line), line);

    return true;
}

LIB_EXPORT rc_t CC SAMExtractorMake(Extractor **state, const KFile * fin, int32_t num_threads=-1)
{
    Extractor * s=(Extractor *)calloc(1,sizeof(*s));
    *state=s;

    VectorInit(&s->headers,0,0);
    VectorInit(&s->alignments,0,0);
    VectorInit(&s->tagvalues,0,0);
    pool_init();

    s->prev_headers=NULL;
    s->prev_aligns=NULL;
    s->pos=0;

    s->hashdvn=false;
    s->hashdso=false;
    s->hashdgo=false;
    s->hassqsn=false;
    s->hassqln=false;
    s->hasrgid=false;
    s->haspgid=false;

    s->rc=0;

    s->infile=fin;

    s->num_threads=num_threads;

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorRelease(Extractor *s)
{
    DBG("release_Extractor");
    SAMlex_destroy();

    SAMExtractorInvalidateHeaders(s);
    SAMExtractorInvalidateAlignments(s);

    VectorWhack(&s->headers,NULL,NULL);

    pool_release();

    memset(s,0,sizeof(Extractor));
    free(s);

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorGetHeaders(Extractor *s, Vector *headers)
{
    rc_t rc=0;
    DBG("get_headers");

    u64 sz=0;
    rc=KFileSize(s->infile,&sz);
    DBG("File size=%'d", sz);
    if (sz < 12)
    {
        ERR("File too small");
        return RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
    }

    rc=readfile(s);
    if (rc) return rc;

    if (!memcmp(readbuf,"\x1f\x8b\x08",3))
    {
        DBG("gzip file, BAM or SAM.gz");
        threadinflate(s);
    } else if (readbuf[0]=='@')
    {
        // TODO: Move to function
        DBG("SAM file");

        while (readline(s))
        {
            rc=SAM_parseline(s);
            if (rc) return rc;

            if (curline[0]!='@')
            {
                // First line of alignments will be processed
                DBG("out of headers");
                break;
            }
        }
        DBG("Done parsing headers");
    }
    else
    {
        ERR("Unkown magic, not a SAM file.");
        return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
    }

    VectorInit(headers,0,0);
    VectorCopy(&s->headers,headers);
    s->prev_headers=headers;
    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorInvalidateHeaders(Extractor *s)
{
    DBG("invalidate_headers");
    for (u32 i=0; i!=VectorLength(&s->headers); ++i)
    {
        Header * hdr=(Header *)VectorGet(&s->headers,i);

        hdr->headercode=NULL;

        Vector * tvs=&hdr->tagvalues;
        for (u32 j=0; j!=VectorLength(tvs); ++j)
        {
            TagValue * tv=(TagValue *)VectorGet(tvs,j);
            tv->tag=NULL;
            tv->value=NULL;
        }
        VectorWhack(&hdr->tagvalues,NULL,NULL);
        hdr=NULL;
    }
    pool_release();
    pool_init();
    VectorWhack(&s->headers,NULL,NULL);
    VectorWhack(&s->tagvalues,NULL,NULL);
    VectorWhack(s->prev_headers,NULL,NULL);
    s->prev_headers=NULL;
    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorGetAlignments(Extractor *s, Vector *alignments)
{
    rc_t rc=0;
    DBG("get_alignments");
    SAMExtractorInvalidateAlignments(s);
    VectorInit(&s->alignments,0,0);

    while (VectorLength(&s->alignments) < 20)
    {
        if (!readline(s)) break;

        if (curline[0]=='@')
        {
            ERR("header restarted");
            break;
        }

        rc=SAM_parseline(s);
        if (rc) return rc;
    }

    DBG("Done parsing %d alignments", VectorLength(&s->alignments));

    VectorInit(alignments,0,0);
    VectorCopy(&s->alignments,alignments);
    VectorWhack(&s->alignments,NULL,NULL);
    s->prev_aligns=alignments;
    DBG("got_alignments");

    return 0;
}

LIB_EXPORT rc_t CC SAMExtractorInvalidateAlignments(Extractor *s)
{
    DBG("invalidate_alignments");
    pool_release();
    pool_init();
    VectorWhack(&s->alignments,NULL,NULL);
    VectorWhack(s->prev_aligns,NULL,NULL);
    s->prev_aligns=NULL;

    return 0;
}

