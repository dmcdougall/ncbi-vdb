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
    #include <strtol.h>
    #include <klib/rc.h>
    #include "samextract.h"
    #include <align/samextract-lib.h>
    int SAMlex(Extractor *);
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

%token TAB
%token <strval> CONTROLCHAR
%token EOL
%token END 0 "end of file"

%expect 0
%%

 /* Bison grammar rules */
sam: /* beginning of input */
   %empty
   | sam line
   ;

line:
    EOL /* Spec is unclear about empty lines, accept for now */
   | CONTROLCHAR { ERR("CONTROLCHAR %d", $1[0]);
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
       COMMENT { mark_headers(state,"CO"); }
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

        mark_headers(state,"HD"); }
    /* TODO: Duplicate header */
    ;

headerlist:   hdr
            | headerlist hdr
  ;

hdr: HDVN VALUE {
        state->hashdvn=true;
        process_header(state,"HD","VN",$2);
        pool_free($2); }
   | HDSO VALUE {
        state->hashdso=true;
        process_header(state,"HD","SO",$2);
        pool_free($2); }
   | HDGO VALUE {
        state->hashdgo=true;
        process_header(state,"HD","GO",$2);
        pool_free($2); }
        /* TODO: Handle >2 tabs in a row */
        /*
  | TAB TAB {
        ERR("two tabs");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
        */
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
        mark_headers(state,"SQ"); }
    ;

sequencelist: sq
    | sequencelist sq
    ;

sq:
      SQSN VALUE {
        state->hassqsn=true;
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | SQLN VALUE {
        if (!inrange($2,1,INT32_MAX))
        {
            ERR("SQ LN field not in range %s",$2);
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        state->hassqln=true;
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | SQAS VALUE {
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | SQM5 VALUE {
        if (!ismd5($2))
            WARN("M5 value not followed by MD5");
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | SQSP VALUE {
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | SQUR VALUE {
        process_header(state,"SQ",$1,$2);
        pool_free($2); }
    | TAB { WARN("Unexpected tab in sequence"); }
    ;

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
        mark_headers(state,"PG"); }
     ;

programlist: pg
    | programlist pg
    ;

pg:
      PGID VALUE {
        state->haspgid=true;
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | PGPN VALUE {
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | PGCL VALUE {
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | PGPP VALUE {
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | PGDS VALUE {
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | PGVN VALUE {
        process_header(state,"PG",$1,$2);
        pool_free($2); }
    | VALUE {
        WARN("Bogus value in PG:%s",$1);
        pool_free($1); }
    ;


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

        mark_headers(state,"RG"); }
    ;

readgrouplist:   rg
            | readgrouplist rg
  ;

rg:  RGID VALUE {
        state->hasrgid=true;
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGCN VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGDS VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGDT VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGFO VALUE {
        if (!isfloworder($2))
            WARN("Flow order incorrec");
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGKS VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGLB VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPG VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPI VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPL PLATFORM {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPL VALUE {
        ERR("Invalid Platform %s", $2);
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPM VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGPU VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | RGSM VALUE {
        process_header(state,"RG",$1,$2);
        pool_free($2); }
   | VALUE VALUE {
        WARN("Unknown readgroup (RG) tag:%s", $1);
        pool_free($1);
        pool_free($2);
        }
   | TAB TAB EOL {
        ERR("empty tags");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
   | TAB EOL { WARN("empty tags"); }
   ;



alignment:
     QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL EOL
     {
        DBG("alignment record");

        u32 iflag = strtou32($2, NULL, 10);
        i32 ipos = strtoi32($4, NULL, 10);

        process_alignment(state,$1,iflag,$3,ipos,$5,$6,$7,$8,$9,$10,$11);

        pool_free($1);
        pool_free($2);
        pool_free($3);
        pool_free($4);
        pool_free($5);
        pool_free($6);
        pool_free($7);
        pool_free($8);
        pool_free($9);
        pool_free($10);
        pool_free($11); }
     | QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL optlist EOL
    {
        DBG("alignment record with optional tags");

        u32 iflag = strtou32($2, NULL, 10);
        i32 ipos = strtoi32($4, NULL, 10);

        process_alignment(state,$1,iflag,$3,ipos,$5,$6,$7,$8,$9,$10,$11);

        pool_free($1);
        pool_free($2);
        pool_free($3);
        pool_free($4);
        pool_free($5);
        pool_free($6);
        pool_free($7);
        pool_free($8);
        pool_free($9);
        pool_free($10);
        pool_free($11); }
    ;

optlist: opt { DBG("opt"); }
       | optlist opt { DBG(" opts"); }
    ;

opt:
    OPTTAG OPTATYPE OPTAVALUE
    {
        DBG("?AA");
        pool_free($1);
        pool_free($3); }
    | OPTTAG OPTITYPE OPTIVALUE
    {
        DBG("?II");
        pool_free($1);
        pool_free($3); }
    | OPTTAG OPTFTYPE OPTFVALUE
    {
        DBG("?FF");
        pool_free($1);
        pool_free($3); }
    | OPTTAG OPTZTYPE OPTZVALUE
    {
        DBG("?ZZ");
        pool_free($1);
        pool_free($3); }
    | OPTTAG OPTHTYPE OPTHVALUE
    {
        DBG("?HH");
        pool_free($1);
        pool_free($3); }
    | OPTTAG OPTBTYPE OPTBVALUE
    {
        DBG("?BB");
        pool_free($1);
        pool_free($3); }
    | OPTITAG OPTITYPE OPTIVALUE
    {
        DBG("III");
        pool_free($1);
        pool_free($3); }
    | OPTZTAG OPTZTYPE OPTZVALUE
    {
        DBG("ZZZ");
        pool_free($1);
        pool_free($3); }
    | OPTBTAG OPTBTYPE OPTBVALUE
    {
        DBG("BBB");
        pool_free($1);
        pool_free($3); }
    ;

%%

