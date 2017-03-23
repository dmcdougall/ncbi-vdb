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

/* %{
   Prologue
   %}
   Declarations
   %%
   Grammar rules
   %%
   Epilogue
   */

%{
    #include <stdio.h>
    #include <ctype.h>
    #include <stdlib.h>
    #include <string.h>
    #include <sys/stat.h>
    #include <sys/types.h>
    #include <regex.h>
    #include <klib/rc.h>
    #include "samextract.h"
    #include <align/samextract-lib.h>
    #include <errno.h>
    #include <strtol.h>

int SAMlex (Extractor *);

void SAMerror(Extractor * state, const char * s)
{
    ERR(" Parsing error: %s",s);
    rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
    state->rc=rc;
    return;
}

/*TODO: Replace with pool allocator. */
static void * myalloc(Extractor * state,size_t sz)
{
    void * buf=malloc(sz);
    if (buf==NULL)
    {
        ERR("out of memory");
        return NULL;
    }
    memset(buf,0,sz);
    VectorAppend(&state->allocs,NULL,buf);
    return buf;
}

static void * mystrdup(Extractor * state,const char * str)
{
    size_t len=strlen(str)+1;
    void * buf=myalloc(state,len);
    memmove(buf,str,len);
    return buf;
}

/* low<=str<=high */
static bool inrange(const char * str, i64 low, i64 high)
{
    i64 i=strtoi64(str, NULL, 10);
    if (errno) return false;
    if (i<low || i>high) return false;
    return true;
}

/* Avoiding handling this as flex token because so many other ReadGroup values
 * could end up looking like flow orders.
 */
static bool isfloworder(const char * str)
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

static rc_t process_header(Extractor * state, const char * type, const char * tag, const char * value)
{
    DBG("processing type %s tag %s value %s", type, tag, value);
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

    TagValue * tv=myalloc(state,sizeof(TagValue));
    if (tv==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    tv->tag=mystrdup(state,tag);
    if (tv->tag==NULL)
    {
        ERR("out of memory");
        rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
        state->rc=rc;
        return rc;
    }
    tv->value=mystrdup(state,value);
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

static rc_t mark_headers(Extractor * state, const char * type)
{
    DBG("mark_headers");
    Header * hdr=(Header *)myalloc(state,sizeof(Header));
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

%}

/* Bison Declarations */
%union {
 char * strval;
}

%name-prefix "SAM"
%param { Extractor * state}
%require "3.0"
%define parse.error verbose

%token <strval> HEADER
%token <strval> SEQUENCE
%token <strval> READGROUP
%token <strval> PROGRAM
%token <strval> COMMENT

%token <strval> VALUE

%token <strval> QNAME
%token <strval> FLAG
%token <strval> RNAME
%token <strval> POS
%token <strval> MAPQ
%token <strval> CIGAR
%token <strval> RNEXT
%token <strval> PNEXT
%token <strval> TLEN
%token <strval> SEQ
%token <strval> QUAL

%token <strval> OPTTAG
%token <strval> OPTITAG
%token <strval> OPTZTAG
%token <strval> OPTBTAG

%token <strval> OPTATYPE
%token <strval> OPTITYPE
%token <strval> OPTFTYPE
%token <strval> OPTZTYPE
%token <strval> OPTHTYPE
%token <strval> OPTBTYPE

%token <strval> OPTAVALUE
%token <strval> OPTIVALUE
%token <strval> OPTFVALUE
%token <strval> OPTZVALUE
%token <strval> OPTHVALUE
%token <strval> OPTBVALUE

%token HDVN
%token HDSO
%token HDGO

%token <strval> RGID
%token <strval> RGCN
%token <strval> RGDS
%token <strval> RGDT
%token <strval> RGFO
%token <strval> RGKS
%token <strval> RGLB
%token <strval> RGPG
%token <strval> RGPI
%token <strval> RGPL
%token <strval> RGPM
%token <strval> RGPU
%token <strval> RGSM
%token <strval> PLATFORM

%token <strval> PGID
%token <strval> PGPN
%token <strval> PGCL
%token <strval> PGPP
%token <strval> PGDS
%token <strval> PGVN

%token <strval> SQSN
%token <strval> SQLN
%token <strval> SQAS
%token <strval> SQM5
%token <strval> SQSP
%token <strval> SQUR

%token COLON
%token TAB
%token CONTROLCHAR
%token EOL
%token END 0 "end of file"

%expect 2
 /* TODO, two shift-reduce conflicts? */
%%

 /* Bison grammar rules */
sam: /* beginning of input */
   %empty
   | sam line
   ;

line:
    EOL /* Spec is unclear about empty lines, accept for now */
   | CONTROLCHAR { ERR("CONTROLCHAR");
                   rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                   state->rc=rc;
   }
   | comment { DBG("comment"); }
   | header { DBG("header done"); }
   | sequence { DBG("sequence"); }
   | program { DBG("program"); }
   | readgroup { DBG("readgroup"); }
   | alignment { DBG("alignment"); }
   ;

comment:
       COMMENT {
       mark_headers(state,"CO");
    }
    ;

header:
      HEADER headerlist EOL
    {
        DBG("header list");
        if (!state->hashdvn)
        {
            ERR("VN tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        if (state->hashdso && state->hashdgo)
           WARN("Both SO and GO tags present");
        if (!state->hashdso && !state->hashdgo)
           WARN("neither SO or GO tags present");

        mark_headers(state,"HD");
        // TODO: Duplicate header
    }
    ;

headerlist:   hdr
            | headerlist hdr
  ;

hdr: HDVN VALUE {
        state->hashdvn=true;
        process_header(state,"HD","VN",$2);
        free($2);
   }
   | HDSO VALUE {
        state->hashdso=true;
        process_header(state,"HD","SO",$2);
        free($2);
   }
   | HDGO VALUE {
        state->hashdgo=true;
        process_header(state,"HD","GO",$2);
        free($2);
   };
  | TAB TAB {
        ERR("two tabs"); /* TODO: Handle >2 tabs in a row */
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
  }
  | TAB { WARN("empty tags"); }
  ;



sequence:
     SEQUENCE sequencelist EOL
    {
        DBG("sequence");
        if (!state->hassqsn)
        {
            ERR("SN tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        if (!state->hassqln)
        {
            ERR("LN tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        mark_headers(state,"SQ");
    }
    ;

sequencelist: sq
    | sequencelist sq
    ;

sq:
      SQSN VALUE {
        state->hassqsn=true;
        process_header(state,"SQ",$1,$2);
        free($2); }
    | SQLN VALUE {
        if (!inrange($2,1,INT32_MAX))
        {
            ERR("SQ LN field not in range %s",$2);
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        state->hassqln=true;
        process_header(state,"SQ",$1,$2);
        free($2); }
    | SQAS VALUE {
        process_header(state,"SQ",$1,$2);
        free($2); }
    | SQM5 VALUE {
        process_header(state,"SQ",$1,$2);
        free($2); }
    | SQSP VALUE {
        process_header(state,"SQ",$1,$2);
        free($2); }
    | SQUR VALUE {
        process_header(state,"SQ",$1,$2);
        free($2); }

program:
      PROGRAM programlist EOL
     {
        DBG("program");
        if (!state->haspgid)
        {
            ERR("ID tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        mark_headers(state,"PG");
     }
     ;

programlist: pg
    | programlist pg
    ;

pg:
      PGID VALUE {
        state->haspgid=true;
        process_header(state,"PG",$1,$2);
        free($2); }
    | PGPN VALUE {
        process_header(state,"PG",$1,$2);
        free($2); }
    | PGCL VALUE {
        process_header(state,"PG",$1,$2);
        free($2); }
    | PGPP VALUE {
        process_header(state,"PG",$1,$2);
        free($2); }
    | PGDS VALUE {
        process_header(state,"PG",$1,$2);
        free($2); }
    | PGVN VALUE {
        process_header(state,"PG",$1,$2);
        free($2); }


readgroup:
      READGROUP readgrouplist EOL
    {
        DBG("readgroup ");
        if (!state->hasrgid)
        {
            ERR("ID tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }

        mark_headers(state,"RG");
    }
    ;

readgrouplist:   rg
            | readgrouplist rg
  ;

rg:  RGID VALUE {
        state->hasrgid=true;
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGCN VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGDS VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGDT VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGFO VALUE {
        if (!isfloworder($2))
            WARN("Flow order incorrec");
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGKS VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGLB VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPG VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPI VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPL PLATFORM {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPL VALUE {
        ERR("Invalid Platform %s", $2);
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPM VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGPU VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | RGSM VALUE {
        process_header(state,"RG",$1,$2);
        free($2); }
   | TAB TAB {
        ERR("two tabs"); /* TODO: Handle >2 tabs in a row */
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
   | TAB TAB EOL {
        ERR("empty tags");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
   | TAB EOL { WARN("empty tags"); }
  ;



alignment:
     QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL optlist EOL
    {
        DBG("Done alignment record");
        Alignment * align=myalloc(state,sizeof(Alignment));
        if (align==NULL)
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            state->rc=rc;
        }
        const char *qname=$1;
        const char *flag=$2;
        const char *rname=$3;
        const char *pos=$4;
        const char *mapq=$5;
        const char *cigar=$6;
        const char *rnext=$7;
        const char *pnext=$8;
        const char *tlen=$9;
        const char *seq=$10;
        const char *qual=$11;
        DBG("align %s %s %s", qname, rnext, qual); // TODO silence warning for now

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

        align->read=mystrdup(state,seq);
        align->cigar=mystrdup(state,cigar);
        align->rname=mystrdup(state,rname);
        align->pos=strtou32(pos,NULL,10);
        align->flags=strtou32(flag,NULL,10);
        VectorAppend(&state->alignments,NULL,align);
        free($1);
        free($2);
        free($3);
        free($4);
        free($5);
        free($6);
        free($7);
        free($8);
        free($9);
        free($10);
        free($11);
    }
    ;

optlist: opt { DBG("opt"); }
       | optlist opt { DBG(" opts"); }
    ;

opt:
    OPTTAG OPTATYPE OPTAVALUE
    {
        DBG("?AA");
        free($1);
        free($3); }
  | OPTTAG OPTITYPE OPTIVALUE
    {
        DBG("?II");
        free($1);
        free($3); }
  | OPTTAG OPTFTYPE OPTFVALUE
    {
        DBG("?FF");
        free($1);
        free($3); }
  | OPTTAG OPTZTYPE OPTZVALUE
    {
        DBG("?ZZ");
        free($1);
        free($3); }
  | OPTTAG OPTHTYPE OPTHVALUE
    {
        DBG("?HH");
        free($1);
        free($3); }
  | OPTTAG OPTBTYPE OPTBVALUE
    {
        DBG("?BB");
        free($1);
        free($3); }
  | OPTITAG OPTITYPE OPTIVALUE
  {
        DBG("III");
        free($1);
        free($3); }
  | OPTZTAG OPTZTYPE OPTZVALUE
  {
        DBG("ZZZ");
        free($1);
        free($3); }
  | OPTBTAG OPTBTYPE OPTBVALUE
  {
        DBG("BBB");
        free($1);
        free($3); }
  ;


%%

