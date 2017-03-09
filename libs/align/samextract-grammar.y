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
    #include <stdint.h>
    #include <klib/rc.h>
    #include "samextract.h"
    #include <align/samextract-lib.h>
    #include "samextract-tokens.h"


size_t alignfields=2; // 1 based, QNAME is #1

int SAMerror(Extractor * state, const char * s)
{
    ERR("Bison error: %s",s);
    rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
    state->rc=rc;
    return rc;
}

//TODO: Replace with pool allocator.
void * myalloc(Extractor * state,size_t sz)
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

void * mystrdup(Extractor * state,const char * str)
{
    size_t len=strlen(str)+1;
    void * buf=myalloc(state,len);
    memmove(buf,str,len);
    return buf;
}

// Returns 1 if match found
int regexcheck(Extractor * state, const char *regex, const char * value)
{
    regex_t preg;

    int result=regcomp(&preg, regex, REG_EXTENDED);
    if (result)
    {
        size_t s=regerror(result, &preg, NULL, 0);
        char *errmsg=malloc(s);
        regerror(result, &preg, errmsg, s);
        ERR("regcomp error on '%s': %s", regex, errmsg);
        free(errmsg);
        regfree(&preg);
        return 0;
    }

    regmatch_t matches[1];
    if (regexec(&preg, value, 1, matches, 0))
    {
        ERR("Value: '%s' doesn't match regex '%s'", value, regex);
        regfree(&preg);
        return 0;
    }
    regfree(&preg);
    return 1;
}

// Returns 1 if OK
int validate(Extractor * state, const char * tag, const char * value)
{
    /* Pair of TAG, regexp: "/..." TODO: or integer range "1-12345" */
    const char * validations[] =
    {
        "VN", "/.*", // @PG also has "/[0-9]+\\.[0-9]+",
        "SO", "/unknown|unsorted|queryname|coordinate",
        "GO", "/none|query|reference",
        "SN", "/[!-)+-<>-~][!-~]*",
        "LN", "/[0]*[1-9][0-9]{0,10}", // TODO: range check 1..2**31-1
        "AS", "/.*",
        "MD", "/[0-9A-Z\\*]{32}", // bam.c treats same as M5
        "M5", "/[0-9A-Za-z\\*]{32}", // TODO: lowercase acceptable?
        "SP", "/.*",
        "UR", "/.*",
        "ID", "/.*",
        "CN", "/.*",
        "DS", "/.*",
        "DT", "/.*",
        "FO", "/\\*|[ACMGRSVTWYHKDBN]+",
        "KS", "/.*",
        "LB", "/.*",
        "PG", "/.*",
        "PI", "/.*",
        "PL", "/.*",
        "PM", "/.*",
        "PU", "/.*",
        "SM", "/.*",
        "PN", "/.*",
        "CL", "/.*",
        "PP", "/.*",
        "DS", "/.*",
        "\0", "\0"
    };

    int ok=0;

    for (size_t i=0;;++i)
        {
            const char *valtag=validations[i*2];
            const char *valval=validations[i*2+1];
            if (*valtag=='\0')
            {
                WARN("No validation for tag %s", tag);
                ok=1;
                break;
            }
            if (!strcmp(tag, valtag))
            {
                if (valval[0]=='/')
                {
                    ok=regexcheck(state,valval+1, value);
                    break;
                } else
                {
                // Parse integer range
                    WARN("range not implemented");
                    ok=1;
                }
            }
        }

    return ok;
}

rc_t check_required_tag(Extractor * state, const char * tags, const char * tag)
{
    if (!strstr(tags,tag))
    {
        ERR("%s tag not seen in header", tag);
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
    }
    return 0;
}

rc_t checkopttagtype(Extractor * state,const char * optfield)
{
    const char *opttypes="AMi ASi BCZ BQZ CCZ CMi COZ CPi CQZ CSZ CTZ E2Z FIi FSZ FZZ H0i H1i H2i HIi IHi LBZ MCZ MDZ MQi NHi NMi OCZ OPi OQZ PGZ PQi PTZ PUZ QTZ Q2Z R2Z RGZ RTZ SAZ SMi TCi U2Z UQi";
    const char type=optfield[3];
    char tag[3];

    tag[0]=optfield[0];
    tag[1]=optfield[1];
    tag[2]='\0';

    if (tag[0]=='X' ||
        tag[0]=='Y' ||
        tag[0]=='Z') return 0;

    const char *p=strstr(opttypes,tag);
    if (p==NULL) return 0;

    if (p[2]!=type)
    {
        ERR("tag %s should have type %c, not %c", tag, p[2], type);
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
    }

    return 0;
}

rc_t process_tagvalue(Extractor * state, const char * tag, const char * value)
{
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
    } else
    {
        if (!validate(state, tag, value))
        {
            ERR("Tag validataion %s failed",tag);
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
            return rc;
        }
        state->tags=realloc(state->tags, strlen(state->tags) + strlen(tag) + 1 + 1);
        strcat(state->tags,tag); strcat(state->tags," ");

        if (!strcmp(tag,"SN"))
        {
            char * s=malloc(strlen(value)+2);
            if (s==NULL)
            {
                ERR("out of memory");
                rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
                }
            strcpy(s,value);
            strcat(s," ");
            if (strstr(state->seqnames,s))
            {
                ERR("duplicate sequence %s", value);
                rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            state->seqnames=realloc(state->seqnames,strlen(state->seqnames) + strlen(value) + 1 + 1);
            if (state->seqnames==NULL)
            {
                ERR("out of memory");
                rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            strcat(state->seqnames,s);
            free(s);
        }
        if (!strcmp(tag,"ID"))
        {
            char * s=malloc(strlen(value)+2);
            if (s==NULL)
            {
                ERR("out of memory");
                rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            strcpy(s,value);
            strcat(s," ");
            if (strstr(state->ids,s))
            {
                ERR("duplicate id %s", value);
                rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            state->ids=realloc(state->ids,strlen(state->ids) + strlen(value) + 1 + 1);
            if (state->ids==NULL)
            {
                ERR("out of memory");
                rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            strcat(state->ids,s);
            free(s);
        }
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
    u32 block=VectorBlock(&state->tagvalues);
    DBG("block is %d",block);
    DBG("Appending %d",VectorLength(&state->tagvalues));

    return 0;
}

rc_t mark_headers(Extractor * state, const char * type)
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

rc_t process_align(Extractor * state, const char *field)
{
    rc_t rc=0;
    const char * opt="(required)";
    if (alignfields>=12) opt="(optional)";
    DBG("rc=%d", state->rc);
    DBG("alignvalue #%zu%s: %s", alignfields, opt, field);
    switch (alignfields)
    {
        case 2: // FLAG
        {
            int flag;
            if (sscanf(field, "%d", &flag)!=1 ||
                flag < 0 ||
                flag > 4095)
            {
                ERR("error parsing FLAG: %s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("flag is %d",flag);
            break;
        }
        case 3: // RNAME
        {
            const char * rname=field;
            if (!regexcheck(state,"\\*|[!-)+-<>-~][!-~]*",rname))
            {
                ERR("error parsing RNAME");
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("rname is %s",rname);
            state->rname=mystrdup(state,rname);
            if (state->rname==NULL)
            {
                ERR("NULL rname");
                rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            break;
        }
        case 4: // POS
        {
            int pos;
            if (sscanf(field, "%d", &pos)!=1 ||
                pos < 0 ||
                pos > INT32_MAX)
            {
                ERR("error parsing POS: %s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("pos is %d",pos);
            state->pos=pos;
            break;
        }
        case 5: // MAPQ
        {
            int mapq;
            if (sscanf(field, "%d", &mapq)!=1 ||
                mapq < 0 ||
                mapq > UINT8_MAX)
            {
                ERR("error parsing MAPQ: %s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("mapq is %d", mapq);
            break;
        }
        case 6: // CIGAR
        {
            const char * cigar=field;
            if (!regexcheck(state,"\\*|([0-9]+[MIDNSHPX=])+",cigar))
            {
                ERR("error parsing cigar");
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("cigar is %s",cigar);
            state->cigar=mystrdup(state,cigar);
            if (state->cigar==NULL)
            {
                ERR("out of memory");
                rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            break;
        }
        case 7: // RNEXT
        {
            const char * rnext=field;
            if (!regexcheck(state,"\\*|=|[!-)+-<>-~][!-~]*",rnext))
            {
                ERR("error parsing rnext");
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("rnext is %s",rnext);
            break;
        }
        case 8: // PNEXT
        {
            int pnext;
            if (sscanf(field, "%d", &pnext)!=1 ||
                pnext < 0 ||
                pnext > INT32_MAX)
            {
                ERR("error parsing PNEXT: %s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("pnext is %d",pnext);
            break;
        }
        case 9: // TLEN
        {
            int tlen;
            if (sscanf(field, "%d", &tlen)!=1 ||
                tlen < INT32_MIN ||
                tlen > INT32_MAX)
            {
                ERR("error parsing TLEN: %s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("tlen is %d", tlen);
            break;
        }
        case 10: // SEQ
        {
            const char * seq=field;
            if (!regexcheck(state,"\\*|[A-Za-z=.]+",seq))
            {
                ERR("error parsing seq");
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("seq is %s",seq);
            state->read=mystrdup(state,seq);
            if (state->read==NULL)
            {
                ERR("out of memory");
                rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                state->rc=rc;
                return rc;
            }
            break;
        }
        case 11: // QUAL
        {
            const char * qual=field;
            if (!regexcheck(state,"[!-~]+",qual))
            {
                ERR("error parsing qual");
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
            }
            DBG("qual is %s", qual);
            DBG("rc=%d", state->rc);
            break;
        }
        default: // Optional
        {
        DBG("optional");
//               /TT:t:
            if ((strlen(field)<5) ||
                field[2]!=':' ||
                field[4]!=':')
                {
                ERR("invald tagtypevalue:%s", field);
                rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                state->rc=rc;
                return rc;
                }
            const char type=field[3];
            if (checkopttagtype(state,field))
            {
                ERR("Optional field tag %s doesn't match type", field);
                WARN("Optional field tag %s doesn't match type", field);
            }
            const char * value=&field[5];
            switch (type)
            {
                case 'A':
                    if (!regexcheck(state,"[!-~]", value))
                    {
                        ERR("value doesn't match A type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                        }
                    break;
                case 'i':
                    if (!regexcheck(state,"[-+]?[0-9]+", value))
                    {
                        ERR("value doesn't match i type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                    }
                    break;
                case 'f':
                    if (!regexcheck(state,"[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", value))
                    {
                        ERR("value doesn't match f type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                    }
                    break;
                case 'Z':
                    if (!regexcheck(state,"[ !-~]*", value))
                    {
                        ERR("value doesn't match Z type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                    }
                    break;
                case 'H':
                    if (!regexcheck(state,"([0-9A-F][0-9A-F])*", value))
                    {
                        ERR("value doesn't match H type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                    }
                    break;
                case 'B':
                    if (!regexcheck(state,"[cCsSiIf](,[-+]?[0-9]*\\.?[0-9]+(eE][-+]?[0-9]+)?)+", value))
                    {
                        ERR("value doesn't match B type:%s",value);
                        rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                        state->rc=rc;
                        return rc;
                    }
                    break;
                default:
                    break;
            }
            DBG("optional field:%s", field);
            break;
        }
    }
    ++alignfields;
    return 0;
}

%}

/* Bison Declarations */
%union {
 int intval;
 char * strval;
 double floatval;
}

%name-prefix "SAM"
/* %define api.pure // was pure-parser */
/* %lex-param   { Extractor * state } */
/* %parse-param { Extractor * state } */
%param { Extractor * state}
%require "3.0"
%define parse.error verbose

%token <strval> HEADER
%token <strval> SEQUENCE
%token <strval> READGROUP
%token <strval> PROGRAM
%token <strval> COMMENT
%token <strval> TAG
%token <strval> VALUE
%token <strval> ALIGNVALUE
%token <strval> QNAME
%token COLON
%token TAB
%token CONTROLCHAR
%token EOL
%token END 0 "end of file"


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
                   return rc;
   }
   | comment { DBG("comment"); }
   | header EOL { DBG("header"); }
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
      HEADER tagvaluelist
    {
        DBG("header tagvaluelist");
        check_required_tag(state,state->tags,"VN");
        if (!strcmp(state->tags,"SO ") &&
            !strcmp(state->tags,"GO "))
           WARN("Both SO and GO tags present");
        if (!(strcmp(state->tags,"SO ") ||
              strcmp(state->tags,"GO ")))
           WARN("neither SO or GO tags present");
        free(state->tags);
        state->tags=strdup("");

        mark_headers(state,"HD");
    }
    ;

sequence:
        SEQUENCE tagvaluelist
    {
        DBG("sequence");
        DBG(" sequences were: %s", state->seqnames);
        check_required_tag(state,state->tags,"SN");
        check_required_tag(state,state->tags,"LN");
        free(state->tags);
        state->tags=strdup("");
        mark_headers(state,"SQ");
    }
    ;

program:
       PROGRAM tagvaluelist
     {
        DBG("ids were: %s", state->ids);
        DBG("program");
        check_required_tag(state,state->tags,"ID");
        free(state->tags);
        state->tags=strdup("");
        mark_headers(state,"PG");
     }
     ;


readgroup:
         READGROUP tagvaluelist
     {
        DBG("readgroup");
        DBG("ids were: %s", state->ids);
        check_required_tag(state,state->tags,"ID");
        free(state->tags);
        state->tags=strdup("");
        mark_headers(state,"RG");
     }
     ;

tagvaluelist: tagvalue { DBG(" one tagvaluelist"); }
            | tagvaluelist tagvalue { DBG(" many tagvaluelist"); }
  ;

tagvalue: TAB TAG COLON VALUE {
        DBG("tagvalue:%s=%s", $2, $4);
        const char * tag=$2;
        const char * value=$4;
        process_tagvalue(state,tag,value);
        free($2);
        free($4);
        };
  | TAB TAB TAG COLON VALUE {
        ERR("two tabs");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
  }
  | TAB TAB EOL {
        ERR("empty tags");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc;
        return rc;
  }
  | TAB TAG TAG {
        const char * tag=$2;
        WARN("malformed TAG:VALUE 'TAB %s(NOT COLON)...'", tag);
        }
  | TAB EOL { WARN("empty tags"); }
  ;

alignment:
         QNAME avlist
    {
        DBG(" avlist qname:%s fields=%zu", $1, alignfields);
        alignfields=2;
        Alignment * align=myalloc(state,sizeof(Alignment));
        if (align==NULL)
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            state->rc=rc;
            return rc;
        }
        align->read=state->read;
        align->cigar=state->cigar;
        align->rname=state->rname;
        align->pos=state->pos;
        VectorAppend(&state->alignments,NULL,align);
        free($1);
    }
    ;

avlist:
      av { DBG(" one av"); }
 |    avlist av { DBG("bison: many avlist"); }
    ;

av:
  TAB ALIGNVALUE
    {
        const char * field=$2;
        rc_t rc=process_align(state,field);
        state->rc=rc;
        free($2);
    }
    ;

%%


/* Epilogue */

