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

#include <kfs/file.h>
#include <klib/rc.h>
#include <klib/defs.h>
#include <klib/vector.h>
#include <klib/text.h>
#include <kproc/queue.h>
#include <kproc/thread.hpp>
#include <kproc/timeout.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strtol.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include "samextract.h"
#include "samextract-tokens.h"
#include <align/samextract-lib.h>

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
          DBG("\t\tParser thread checking parser queue");
          if (c && c->state!=uncompressed)
          {
              ERR("\t\tparser queue bad state");
              return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
          }
          if (c) c->state=empty;
          c=NULL;

          while (true)
          {
              TimeoutInit(&tm, 5000); // 5 seconds
              rc_t rc=KQueuePop(que, &where, &tm);
              if (rc==0)
              {
                  c=(chunk *)where;
                  if (c->state==empty)
                  {
                      ERR("\t\tParser queue bad state");
                      return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
                  }
                  while (c->state!=uncompressed)
                  {
                      DBG("\t\tParser busy");
                      usleep(1);
                  }
                  DBG("\t\tParser thread chunk %p size %u", c->in, c->insize);
                  break;
              } else if ((int)GetRCObject(rc)== rcTimeout)
              {
                  INFO("\t\tParser thread queue empty");
                  if (KQueueSealed(que))
                  {
                      INFO("\t\tQueue sealed, Parser thread complete");
                      return false;
                  }
              } else if ((int)GetRCObject(rc)== rcData)
              {
                  DBG("\t\tParser thread queue data (empty?)");
                  return false;
              } else
              {
                  ERR("Parser rc=%d",rc);
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


static rc_t inflater(const KThread * kt, void * in)
{
    Extractor * state=(Extractor *)in;
    struct timeout_t tm;

    z_stream strm;
    pthread_t threadid=pthread_self();
    DBG("\tInflater thread %lu started.",threadid);

    while (true)
    {
        void * where=NULL;
        DBG("\t\tthread %lu checking queue",threadid);
        TimeoutInit(&tm, 5000); // 5 seconds
        rc_t rc=KQueuePop(state->inflatequeue, &where, &tm);
        if (rc==0)
        {
            chunk * c=(chunk *)where;
            DBG("\t\tinflater thread %lu chunk %p size %u", threadid, c->in, c->insize);
            if (c->state!=compressed)
            {
                ERR("Inflater bad state");
                return RC(rcAlign,rcFile,rcReading,rcData,rcInvalid);
            }

            memset(&strm,0,sizeof strm);
            DBG("\tinflating %d bytes", c->insize);
            strm.next_in=c->in;
            strm.avail_in=c->insize;
            strm.next_out=c->out;
            strm.avail_out=c->outsize;
            int zrc=inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
            switch (zrc)
            {
              case Z_OK:
                  break;
              case Z_MEM_ERROR:
                  ERR("Out of memory in zlib");
                  return RC(rcAlign,rcFile,rcReading,rcMemory,rcExhausted);
              case Z_VERSION_ERROR:
                  ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
              case Z_STREAM_ERROR:
                  ERR("zlib stream error");
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
              default:
                  ERR("zlib error %s",strm.msg);
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
            }

            zrc=inflate(&strm,Z_FINISH);
            switch (zrc)
            {
              case Z_OK:
              case Z_STREAM_END:
                  DBG("\t\tthread %lu OK %d %d %lu", threadid, strm.avail_in, strm.avail_out,strm.total_out);
                  c->outsize=strm.total_out;
                  c->state=uncompressed;
                  while (true)
                  {
                    struct timeout_t tm;
                    TimeoutInit(&tm, 1000); // 1 second
                    rc_t rc=KQueuePush(state->parsequeue, (void *)c, &tm);
                    if ((int)GetRCObject(rc)== rcTimeout)
                    {
                        DBG("parse queue full");
                    } else if (rc == 0)
                    {
                        DBG("parse queued: %p %d %d %d", c->in, c->insize, rc, rcTimeout);
                       break;
                    } else
                    {
                        DBG("parse queue %d", rc);
                    }
                  }
                  break;
              case Z_MEM_ERROR:
                  ERR("error: Out of memory in zlib");
                  return RC(rcAlign,rcFile,rcReading,rcMemory,rcExhausted);
              case Z_VERSION_ERROR:
                  ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
              case Z_STREAM_ERROR:
                  ERR("zlib stream error %s",strm.msg);
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
              default:
                  ERR("inflate error %d %s",zrc, strm.msg);
                  return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
            }
            inflateEnd(&strm);
        } else if ((int)GetRCObject(rc)== rcTimeout)
        {
            DBG("\t\tthread %lu queue empty",threadid);
            if (KQueueSealed(state->inflatequeue))
            {
                DBG("\t\tQueue sealed, thread %lu complete", threadid);
                INFO("\t\tinflater thread %lu terminating.",threadid);
                KQueueSeal(state->parsequeue);
                return 0;
            }
        } else if ((int)GetRCObject(rc)== rcData)
        {
            DBG("\t\tthread %lu queue data (empty?)",threadid);
            break;
        } else
        {
            ERR("rc=%d",rc);
            return rc;
        }
    }

    DBG("\t\tinflater thread %lu wrongly terminating.",threadid);
    KQueueSeal(state->parsequeue);
    return 0;
}


static rc_t parser(const KThread * kt, void * in)
{
    Extractor * state=(Extractor *)in;

    pthread_t threadid=pthread_self();
    DBG("\tParser  thread %lu started.",threadid);

    chunkview cv(state->parsequeue);

    char magic[4];
    if (!cv.getbytes(magic,4)) return 1;
    if (memcmp(magic,"BAM\x01",4))
    {
        ERR("BAM magic not found");
        return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
    }
    i32 l_text;
    if (!cv.getbytes(&l_text,4)) return 1;
    DBG("l_text=%d",l_text);
    if (l_text<0)
    {
        ERR("error: invalid l_text");
        return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
    }

    char *text=(char *)calloc(1,l_text + 2);
    if (!cv.getbytes(text,l_text)) return 1;
    text[l_text+1]='\0';

    DBG("SAM header:'%s'",text);
    char * t=text;
    while (strlen(t))
    {
        char * nl=(char *)strchr(t, '\n');
        if (!nl)
        {
            size_t linelen=strlen(t);
            DBG("noln linelen %d",linelen);
            memmove(curline,t,linelen);
            curline[linelen+1]='\n';
            curline[linelen+2]='\0';
            t+=linelen;
        }
        else
        {
            size_t linelen=1+nl-t;
            DBG("ln   linelen %d",linelen);
            memmove(curline,t,linelen);
            curline[linelen]='\0';
            t+=linelen;
        }
        DBG("curline is '%s'",curline);
        if (curline[0]!='@')
            ERR("Not a SAM header line");
        curline_len=strlen(curline);
        SAMparse(state);
    }
    free(text);
    text=NULL;

    DBG("BAM done with headers");
    state->file_status=headers;
    while (state->file_status!=alignments)
    {
        DBG("Waiting for alignments");
        sleep(1);
    }

    i32 n_ref;
    if (!cv.getbytes(&n_ref,4)) return 1;
    if (n_ref<0)
    {
        ERR("error: invalid n_ref");
        return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
    }
    DBG("# references %d", n_ref);

    for (int i=0; i!=n_ref; ++i)
    {
        i32 l_name;
        if (!cv.getbytes(&l_name,4)) return 1;
        DBG("%d: l_name=%d",i, l_name);
        if (l_name < 0)
        {
            ERR("error: invalid reference name length");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }
        if (l_name > 256)
        {
            WARN("warning: Long reference name");
            return 0;
        }
        char *name=(char *)calloc(1,l_name+1);
        if (!cv.getbytes(name,l_name)) return 1;
        DBG("%d: reference name %s",i, name);
        free(name); // TODO, persist?
        name=NULL;
        i32 l_ref;
        if (!cv.getbytes(&l_ref,4)) return 1;
        DBG("length of reference sequence %d=%d",i,l_ref);
    }

    bamalign align;
    while (cv.getbytes(&align,sizeof(align)))
    {
        DBG("alignment block_size=%d refID=%d pos=%d",
            align.block_size,
            align.refID,
            align.pos);

        if (align.block_size < 0)
        {
            ERR("error: invalid block_size");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }
        if (align.pos < 0)
        {
            ERR("error: invalid pos");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }
        if (align.refID < -1 || align.refID > n_ref)
        {
            ERR("error: bad refID");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }
        if (align.next_refID < -1 || align.next_refID > n_ref)
        {
            ERR("error: bad next_refID");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }
        DBG("align.bin_mq_nl=%d",align.bin_mq_nl);
        u16 bin=align.bin_mq_nl >> 16;
        u8 mapq=(align.bin_mq_nl >> 8) & 0xff;
        u8 l_read_name=align.bin_mq_nl & 0xff;
        DBG("bin=%d mapq=%d l_read_name=%d", bin, mapq, l_read_name);

        u16 flag=align.flag_nc >> 16;
        u16 n_cigar_op=align.flag_nc & 0xffff;
        DBG("flag=%x n_cigar_op=%d", flag, n_cigar_op);

        char * read_name=(char *)calloc(1,l_read_name);
        if (!cv.getbytes(read_name,l_read_name)) return 1;
        DBG("read_name=%s",read_name);

        static const char opmap[]="MIDNSHP=X???????";
        u32 * cigar=(u32 *)calloc(n_cigar_op,sizeof(u32));

        if (!cv.getbytes(cigar,n_cigar_op*4)) return 1;
        // Compute size of expanded CIGAR
        size_t rleopslen=0;
        for (int i=0; i!=n_cigar_op; ++i)
        {
            i32 oplen=cigar[i] >> 4;
            if (!oplen) ERR("Bogus CIGAR op length");
            rleopslen+=oplen;
        }
        char * scigar=(char *)calloc(rleopslen+1,1);
        char * p=scigar;
        for (int i=0; i!=n_cigar_op; ++i)
        {
            i32 oplen=cigar[i] >> 4;
            i32 op=cigar[i] & 0xf;
            DBG("\tcigar %d=%x len=%d %d(%c)", i, cigar[i], oplen, op, opmap[op]);
            for (int j=0; j!=oplen; ++j)
                *p++=(char)opmap[op];
        }
        DBG("scigar is '%s'", scigar);

        static const char seqmap[]="=ACMGRSVTWYHKDBN";
        char * seq=(char *)calloc(1,align.l_seq+1);
        int bytesofseq=(align.l_seq+1)/2;
        u8 * seqbytes=(u8 *)calloc(1,bytesofseq);
        if (!cv.getbytes(seqbytes,bytesofseq)) return 1;
        int i=0;
        int j=0;
        while (i < align.l_seq)
        {
            seq[i]=seqmap[seqbytes[j] >> 4];
            i+=2;
            j+=1;
        }
        i=1;
        j=0;
        while (i < align.l_seq)
        {
            seq[i]=seqmap[seqbytes[j] & 0xf];
            i+=2;
            j+=1;
        }

        char *qual=(char*)calloc(1,align.l_seq);
        if (!cv.getbytes(qual,align.l_seq)) return 1;

        DBG("%d pairs in sequence",align.l_seq);
//        for (int i=0; i!=align.l_seq; ++i)
//            DBG(" seq#%d %c %.2x ",i, seq[i], qual[i]);
        DBG("seq=%s",seq);

        int remain=align.block_size-(sizeof(align)+l_read_name+n_cigar_op*4+bytesofseq+align.l_seq)+4; // TODO, why 4?
        DBG("%d bytes remaining for ttvs",remain);
        char * ttvs=(char*)calloc(1,remain);
        if (!cv.getbytes(ttvs,remain)) return 1;
        char * cur=ttvs;
        while (cur<ttvs+remain)
        {
            char tag[2];
            char c;
            i8 i8;
            u8 u8;
            i16 i16;
            u16 u16;
            i32 i32;
            u32 u32;
            char * z;
            tag[0]=*cur++;
            tag[1]=*cur++;
            char val_type=*cur++;
            DBG("ttv: %c%c:%c", tag[0], tag[1], val_type);
            switch (val_type)
            {
              case 'A':
                  c=*cur++;
                  DBG("val='%c'",c);
                  break;
              case 'c':
                  i8=*cur++;
                  DBG("val=%d",i8);
                  break;
              case 'C':
                  u8=*cur++;
                  DBG("val=%d",u8);
                  break;
              case 's':
                  memmove(&i16,cur,2);
                  DBG("val=%d",i16);
                  cur+=2;
                  break;
              case 'S':
                  memmove(&u16,cur,2);
                  DBG("val=%d",u16);
                  cur+=2;
                  break;
              case 'i':
                  memmove(&i32,cur,4);
                  cur+=4;
                  break;
              case 'I':
                  memmove(&u32,cur,4);
                  cur+=4;
                  break;
              case 'f':
                  //float f;
                  break;
              case 'Z':
                  z=cur;
                  while (isprint(*cur))
                      ++cur;
                  DBG("val='%s'",z);
                  ++cur;
                  break;
              case 'H':
                  z=cur;
                  while (isalnum(*cur))
                      ++cur;
                  DBG("val='%s'",z);
                  // TODO: Convert to ?
                  ++cur;
                  break;
              case 'B':
                  val_type=*cur++;
                  memmove(&u32,cur,4);
                  cur+=4;
                  cur+=u32*1; // TODO, based on size of val_type
                  break;
              default:
                  ERR("Bad val_type:%c", val_type);
                  return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
            }
        }
        DBG("no more ttvs");

        // We want read (seq), cigar, rname, pos and flags
        // name=Qname, reference sequence name
        // read_name=qname, query template name

        char sflag[16];
        snprintf(sflag, sizeof sflag, "%u", flag);
        char spos[16];
        snprintf(spos, sizeof spos, "%u", align.pos);

        DBG("sflag='%s', spos='%s'",sflag,spos);
        process_alignment(state,
                          NULL, // QNAME
                          sflag,
                          read_name,  // RNAME
                          spos,
                          "0", //mapq
                          scigar, // cigar
                          NULL, //rnext
                          "0", //pnext
                          "0", //tlen
                          seq, // read
                          NULL //qual
                          );

        free(read_name);
        read_name=NULL;
        free(cigar);
        cigar=NULL;
        free(seq);
        seq=NULL;
        free(seqbytes);
        seqbytes=NULL;
        free(qual);
        qual=NULL;
        free(ttvs);
        ttvs=NULL;
        free(scigar);
        scigar=NULL;
    }
    if (!KQueueSealed(state->parsequeue))
    {
        ERR("out of data but queue not sealed");
        return 1;
    }
    INFO("parser thread complete");
    return 0;
} // parser

static rc_t seeker(const KThread * kt, void * in)
{
    Extractor * state=(Extractor *)in;

    pthread_t threadid=pthread_self();
    DBG("\tSeeker thread %lu started.",threadid);

    state->file_pos=0;
    while (true)
    {
        if (state->readbuf_sz < 28)
        {
            ERR("Small block");
        }

        if (!memcmp(state->readbuf, "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00", 28))
        {
            INFO("EOF marker found");
            break;
        }

        z_stream strm;
        memset(&strm,0,sizeof strm);
        strm.next_in=(Bytef*)state->readbuf;
        strm.avail_in=(uInt)state->readbuf_sz;
        int zrc=inflateInit2(&strm, MAX_WBITS + 16); // Only gzip format
        switch (zrc)
        {
          case Z_OK:
              break;
          case Z_MEM_ERROR:
              ERR("error: Out of memory in zlib");
              return RC(rcAlign,rcFile,rcReading,rcMemory,rcExhausted);
          case Z_VERSION_ERROR:
              ERR("zlib version is not compatible; need version %s but have %s", ZLIB_VERSION,zlibVersion());
              return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
          case Z_STREAM_ERROR:
              ERR("zlib stream error");
              return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
          default:
              ERR("zlib error");
              return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
        }

        gz_header head;
        u8 extra[256];
        memset(&head,0,sizeof head);
        head.extra=extra;
        head.extra_max=sizeof(extra);
        char outbuf[64];
        strm.next_out=(Bytef*)outbuf;
        strm.avail_out=64;
        zrc=inflateGetHeader(&strm,&head);
        while (head.done==0)
        {
            DBG("inflating gzip header");
            zrc=inflate(&strm,Z_BLOCK);
            if (zrc!=Z_OK)
            {
                for (int i=0; i!=4; ++i)
                    DBG("readbuf: %x", (unsigned char)state->readbuf[i]);
                ERR("inflate error %d %s", zrc, strm.msg);
                return RC(rcAlign,rcFile,rcConstructing,rcNoObj,rcUnexpected);
            }
        }
        DBG("found gzip header");
        // BC 02 bb
        if (head.extra && head.extra_len==6 && head.extra[0]=='B' && head.extra[1]=='C' && head.extra[2]==2 && head.extra[3]==0)
        {
            u16 bsize=head.extra[4]+head.extra[5]*256;
            inflateEnd(&strm);
            DBG("total_in:%d",strm.avail_in);
            if (bsize<=28)
            {
                DBG("small block found");
                break;
            }

            size_t block_size=12; // Up to and including XLEN
            block_size+=head.extra_len;
            block_size+=(bsize-head.extra_len-19); // CDATA
            block_size+=8; // CRC32 and isize
            DBG("block_size is %d bsize is %d", block_size, bsize);

            state->file_pos+=block_size;

            chunk * c=(chunk *)calloc(1,sizeof(chunk));
            c->insize=bsize+1;
            memmove(c->in, state->readbuf, block_size+1);
            c->outsize=sizeof(c->out);
            c->state=compressed;

//            VectorAppend(&state->chunks,NULL,c);
            while (true)
            {
                struct timeout_t tm;
                TimeoutInit(&tm, 1000); // 1 second
                rc_t rc=KQueuePush(state->inflatequeue, (void *)c, &tm);
                if ((int)GetRCObject(rc)== rcTimeout)
                {
                    DBG("inflate queue full");
                } else if (rc == 0)
                {
                    DBG("inflate queued: %p %d %d %d", c->in, c->insize, rc, rcTimeout);
                    break;
                } else
                {
                    DBG("inflate queue %d", rc);
                }
            }
        } else
        {
            ERR("error: BAM required extra extension not found");
            return RC(rcAlign,rcFile,rcParsing,rcData,rcInvalid);
        }

        DBG("Getting more");

        DBG("reading in at %d",state->file_pos);
        rc_t rc=KFileReadAll(state->infile, state->file_pos, state->readbuf, state->readbuf_sz, &state->readbuf_sz);
        //state->file_pos+=state->readbuf_sz;
        if (rc) { ERR("readfile error"); return rc;}
        DBG("Read in %d",state->readbuf_sz);
        state->readbuf_pos=0;
        if (!state->readbuf_sz)
        {
            DBG("Buffer complete. EOF");
        }
    }

    KQueueSeal(state->inflatequeue);
    INFO("seeker thread complete");
    return 0;
}

void waitforthreads(Vector * threads)
{
    for (u32 i=0; i!=VectorLength(threads); ++i)
    {
        rc_t rc=0;
        DBG("waiting for thread %d to complete", i);
        KThreadWait((KThread*)VectorGet(threads,i),&rc);
    }
    DBG("all threads completed");
}

void releasethreads(Vector * threads)
{
    DBG("Release threads");
    for (u32 i=0; i!=VectorLength(threads); ++i)
    {
        KThreadRelease((KThread*)VectorGet(threads,i));
    }
}

rc_t threadinflate(Extractor * state)
{
    // Seeker thread
    KThread * seekerthread=(KThread *)calloc(1,sizeof(KThread));
    KThreadMake(&seekerthread, seeker, (void*)state);
    VectorAppend(&state->threads, NULL, seekerthread);

    sleep(1);

    // Inflater thread
    KThread * inflaterthread=(KThread *)calloc(1,sizeof(KThread));
    KThreadMake(&inflaterthread, inflater, (void *)state);
    VectorAppend(&state->threads, NULL, inflaterthread);

    sleep(1);

    // Parse thread
    KThread * parserthread=(KThread *)calloc(1,sizeof(KThread));
    KThreadMake(&parserthread, parser, (void *)state);
    VectorAppend(&state->threads, NULL, parserthread);

    return 0;
}

