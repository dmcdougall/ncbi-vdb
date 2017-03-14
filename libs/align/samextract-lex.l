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

   /*
      lex file is three sections:
      definitions: name definition
      %%
      rules: pattern action
      %%
      user code
   */

  /*
    @HD\tVN:1.4\tSO:coordinate
    @SQ\tSN:1\t...
    @SQ\tSN:16\t...
    @RG\tID:PM89\tPL:Illumina\tLB...
    HWI-...
  */

%top{
#include <stdio.h>
#include <string.h>
#include "samextract.h"
#include "samextract-tokens.h"

#define LEXDEBUG 0

char * lexbuf=NULL;
size_t lexbuf_sz=0;

int moredata(char * buf, int * numbytes, size_t maxbytes);
#define YY_INPUT(buf,numbytes,maxbytes) moredata(buf,&numbytes,maxbytes)


  }

 /* definitions */
digits [[:digit:]]+
integer [-+]?[0-9]+
float  [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
eol    (\r\n|\r|\n)
controlchar [\x00-\x08\x0b-\x0c\x0e-\x1f]

tag [A-Za-z][A-Za-z0-9]

 /* No @ in first, to distinguish from headers */
qname [!-?A-~][!-~]{0,253}

 /* See bison's lex-param/parse-param */
%option prefix="SAM"
 /* stop after first EOF */
%option noyywrap
%option nodefault
%option warn
%option nodebug
 /* Don't use these functions */
%option nounput
%option noinput
 /* Not a TTY */
%option never-interactive
 /* %option yylineno */
%option nounistd

 /* Start states */
%x INALIGNMENT
%x AFTERCOLON

%%
 /* rules */
^@HD   { DBG( " Header"); return HEADER; }
^@SQ   { DBG( " SQ"); return SEQUENCE; }
^@RG   { DBG( " RG"); return READGROUP; }
^@PG   { DBG( " PG"); return PROGRAM; }
^@CO\t.*{eol} { DBG( " Comment"); return COMMENT; }


<INITIAL,INALIGNMENT>^{qname} { BEGIN INALIGNMENT; SAMlval.strval=strdup(yytext); DBG(" alignment qname, INALIGNMENT"); return QNAME; }
<INITIAL,INALIGNMENT,AFTERCOLON>{controlchar} { DBG("CONTROLCHAR1"); return CONTROLCHAR; }
<INITIAL,INALIGNMENT>\t { DBG(" tab"); return TAB;}
<AFTERCOLON>\t { BEGIN INITIAL; DBG(" tab INITIAL"); return TAB;}

<INITIAL,INALIGNMENT,AFTERCOLON>{eol} { BEGIN INITIAL; DBG(" eol INITIAL"); return EOL;}
<INALIGNMENT>[ -~]+ {
    SAMlval.strval=strdup(yytext); DBG(" alignvalue"); return ALIGNVALUE; }
{tag}  {
    SAMlval.strval=strdup(yytext); DBG(" Valid Tag:%s", yytext); return TAG; }
<INITIAL,INALIGNMENT>: { BEGIN AFTERCOLON; DBG(" colon AFTERCOLON"); return COLON;}
<AFTERCOLON>[ -~]+ { SAMlval.strval=strdup(yytext); DBG(" Value:%s", yytext); return VALUE; }

.|\n       DBG("DEFAULT '%c' ", *yytext);
 /* TODO: <<EOF>> rule */
%%

#if YY_FLEX_MAJOR_VERSION == 2 && YY_FLEX_MINOR_VERSION < 6
#error "Flex version should be >= 2.6.X"
#endif


rc_t SAM_parsebegin(Extractor * state)
{
    DBG("parsebegin");
    return 0;
}

rc_t SAM_parseend(Extractor * state)
{
    DBG("parseend");
    free(lexbuf);
    lexbuf=NULL;
    lexbuf_sz=0;
    SAMlex_destroy();
    return 0;
}

rc_t SAM_parsebuffer(Extractor * state, char * str, size_t size)
{
    DBG("Parsingbuffer");
    state->rc=0;
    if (lexbuf==NULL || size > lexbuf_sz)
    {
        lexbuf=realloc(lexbuf,size+2);
        lexbuf_sz=size+2;
        if (lexbuf==NULL)
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign,rcFile,rcReading,rcMemory,rcExhausted);
            state->rc=rc;
            return rc;
        }
    }
    memmove(lexbuf,str,size);
    lexbuf[size]='\0';
    SAMparse(state);
    return state->rc;
}

int moredata(char * buf, int * numbytes, size_t maxbytes)
{
    DBG("moredata %p %d", buf, maxbytes);
    strcpy(buf,lexbuf);
    *numbytes=strlen(buf);
    lexbuf[0]='\0';
    return 0;
}


