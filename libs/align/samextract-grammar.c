/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         SAMparse
#define yylex           SAMlex
#define yyerror         SAMerror
#define yydebug         SAMdebug
#define yynerrs         SAMnerrs

#define yylval          SAMlval
#define yychar          SAMchar

/* Copy the first part of user declarations.  */
#line 37 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:339  */

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

    #define YYDEBUG 1
/*    #define SAMdebug 1 */

    size_t alignfields=2; // 1 based, QNAME is #1

    int SAMerror(Extractor * extractor, const char * s)
    {
        ERR("Bison error: %s",s);
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        globstate->rc=rc;
        extractor->rc=rc;
        return rc;
    }

    void * myalloc(size_t sz)
    {
        void * buf=malloc(sz);
        if (buf==NULL) 
        {
            ERR("out of memory");
            return NULL;
        }
        memset(buf,0,sz);
        VectorAppend(&globstate->allocs,NULL,buf);
        return buf;
    }

    void * mystrdup(const char * str)
    {
        size_t len=strlen(str)+1;
        void * buf=myalloc(len);
        memmove(buf,str,len);
        return buf;
    }
/*
    void * myrealloc(void * ptr, size_t sz)
    {
        for (u32 i=0; i!=VectorLength(&globstate->allocs); ++i)
        {
            void * p=VectorGet(&globstate->allocs,i);
            if (p==ptr) ...

    }
*/
    // Returns 1 if match found
    int regexcheck(const char *regex, const char * value)
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
    int validate(const char * tag, const char * value)
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
                    ok=regexcheck(valval+1, value);
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

    rc_t check_required_tag(const char * tags, const char * tag)
    {
        if (!strstr(tags,tag))
        {
            ERR("%s tag not seen in header", tag);
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            globstate->rc=rc;
            return rc;
        }
        return 0;
    }

    rc_t checkopttagtype(const char * optfield)
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
            globstate->rc=rc;
            return rc;
        }

        return 0;
    }

    rc_t process_tagvalue(const char * tag, const char * value)
    {
        if (strlen(tag)!=2)
        {
            ERR("tag '%s' must be 2 characters", tag);
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            globstate->rc=rc;
            return rc;
        }

        if (islower(tag[0] &&
            islower(tag[1])))
        {
            DBG("optional tag");
        } else
        {
            if (!validate(tag, value))
            {
                ERR("Tag validataion %s failed",tag);
                rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                globstate->rc=rc;
                return rc;
            }
            globstate->tags=realloc(globstate->tags, strlen(globstate->tags) + strlen(tag) + 1 + 1);
            strcat(globstate->tags,tag); strcat(globstate->tags," ");

            if (!strcmp(tag,"SN"))
            {
                char * s=malloc(strlen(value)+2);
                if (s==NULL) 
                {
                    ERR("out of memory");
                    rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                    }
                strcpy(s,value);
                strcat(s," ");
                if (strstr(globstate->seqnames,s))
                {
                    ERR("duplicate sequence %s", value);
                    rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                globstate->seqnames=realloc(globstate->seqnames,strlen(globstate->seqnames) + strlen(value) + 1 + 1);
                if (globstate->seqnames==NULL) 
                {
                    ERR("out of memory");
                    rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                }
                strcat(globstate->seqnames,s);
                free(s);
            }
            if (!strcmp(tag,"ID"))
            {
                char * s=malloc(strlen(value)+2);
                if (s==NULL) 
                {
                    ERR("out of memory");
                    rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                }
                strcpy(s,value);
                strcat(s," ");
                if (strstr(globstate->ids,s))
                {
                    ERR("duplicate id %s", value);
                    rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                globstate->ids=realloc(globstate->ids,strlen(globstate->ids) + strlen(value) + 1 + 1);
                if (globstate->ids==NULL) 
                {
                    ERR("out of memory");
                    rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                }
                strcat(globstate->ids,s);
                free(s);
            }
        }

        TagValue * tv=myalloc(sizeof(TagValue));
        if (tv==NULL) 
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            globstate->rc=rc;
            return rc;
        }
        tv->tag=mystrdup(tag);
        if (tv->tag==NULL) 
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            globstate->rc=rc;
            return rc;
        }
        tv->value=mystrdup(value);
        if (tv->value==NULL) 
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            globstate->rc=rc;
            return rc;
        }
        VectorAppend(&globstate->tagvalues,NULL,tv);
        u32 block=VectorBlock(&globstate->tagvalues);
        DBG("block is %d",block);
        DBG("Appending %d",VectorLength(&globstate->tagvalues));
        return 0;
    }

    rc_t mark_headers(const char * type)
    {
        DBG("mark_headers");
        Header * hdr=(Header *)myalloc(sizeof(Header));
        if (hdr==NULL) 
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            globstate->rc=rc;
            return rc;
        }
        hdr->headercode=type;
        VectorCopy(&globstate->tagvalues,&hdr->tagvalues);
        VectorAppend(&globstate->headers,NULL,hdr);
        VectorWhack(&globstate->tagvalues,NULL,NULL);
        return 0;
    }

    rc_t process_align(const char *field)
    {
        rc_t rc=0;
        const char * opt="(required)";
        if (alignfields>=12) opt="(optional)";
        DBG("rc=%d", globstate->rc);
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
                    globstate->rc=rc;
                    return rc;
                }
                DBG("flag is %d",flag);
                break;
            }
            case 3: // RNAME
            {
                const char * rname=field;
                if (!regexcheck("\\*|[!-)+-<>-~][!-~]*",rname))
                {
                    ERR("error parsing RNAME");
                    rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                DBG("rname is %s",rname);
                globstate->rname=mystrdup(rname);
                if (globstate->rname==NULL) 
                {
                    ERR("NULL rname");
                    rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
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
                    globstate->rc=rc;
                    return rc;
                }
                DBG("pos is %d",pos);
                globstate->pos=pos;
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
                    globstate->rc=rc;
                    return rc;
                }
                DBG("mapq is %d", mapq);
                break;
            }
            case 6: // CIGAR
            {
                const char * cigar=field;
                if (!regexcheck("\\*|([0-9]+[MIDNSHPX=])+",cigar))
                {
                    ERR("error parsing cigar");
                    rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                DBG("cigar is %s",cigar);
                globstate->cigar=mystrdup(cigar);
                if (globstate->cigar==NULL) 
                {
                    ERR("out of memory");
                    rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                }
                break;
            }
            case 7: // RNEXT
            {
                const char * rnext=field;
                if (!regexcheck("\\*|=|[!-)+-<>-~][!-~]*",rnext))
                {
                    ERR("error parsing rnext");
                    rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
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
                    globstate->rc=rc;
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
                    globstate->rc=rc;
                    return rc;
                }
                DBG("tlen is %d", tlen);
                break;
            }
            case 10: // SEQ
            {
                const char * seq=field;
                if (!regexcheck("\\*|[A-Za-z=.]+",seq))
                {
                    ERR("error parsing seq");
                    rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                DBG("seq is %s",seq);
                globstate->read=mystrdup(seq);
                if (globstate->read==NULL) 
                {
                    ERR("out of memory");
                    rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
                    globstate->rc=rc;
                    return rc;
                }
                break;
            }
            case 11: // QUAL
            {
                const char * qual=field;
                if (!regexcheck("[!-~]+",qual))
                {
                    ERR("error parsing qual");
                    rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                    globstate->rc=rc;
                    return rc;
                }
                DBG("qual is %s", qual);
                DBG("rc=%d", globstate->rc);
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
                    globstate->rc=rc;
                    return rc;
                  }
                const char type=field[3];
                if (checkopttagtype(field))
                {
                    ERR("Optional field tag %s doesn't match type", field);
                    WARN("Optional field tag %s doesn't match type", field);
                }
                const char * value=&field[5];
                switch (type)
                {
                    case 'A':
                        if (!regexcheck("[!-~]", value))
                        {
                            ERR("value doesn't match A type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
                            return rc;
                            }
                        break;
                    case 'i':
                        if (!regexcheck("[-+]?[0-9]+", value))
                        {
                            ERR("value doesn't match i type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
                            return rc;
                        }
                        break;
                    case 'f':
                        if (!regexcheck("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", value))
                        {
                            ERR("value doesn't match f type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
                            return rc;
                        }
                        break;
                    case 'Z':
                        if (!regexcheck("[ !-~]*", value))
                        {
                            ERR("value doesn't match Z type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
                            return rc;
                        }
                        break;
                    case 'H':
                        if (!regexcheck("([0-9A-F][0-9A-F])*", value))
                        {
                            ERR("value doesn't match H type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
                            return rc;
                        }
                        break;
                    case 'B':
                        if (!regexcheck("[cCsSiIf](,[-+]?[0-9]*\\.?[0-9]+(eE][-+]?[0-9]+)?)+", value))
                        {
                            ERR("value doesn't match B type:%s",value);
                            rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                            globstate->rc=rc;
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


#line 668 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "samextract-tokens.h".  */
#ifndef YY_SAM_HOME_VARTANIANMH_DEVEL_NCBI_VDB_LIBS_ALIGN_SAMEXTRACT_TOKENS_H_INCLUDED
# define YY_SAM_HOME_VARTANIANMH_DEVEL_NCBI_VDB_LIBS_ALIGN_SAMEXTRACT_TOKENS_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int SAMdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    END = 0,
    HEADER = 258,
    SEQUENCE = 259,
    READGROUP = 260,
    PROGRAM = 261,
    COMMENT = 262,
    TAG = 263,
    VALUE = 264,
    ALIGNVALUE = 265,
    QNAME = 266,
    COLON = 267,
    TAB = 268,
    CONTROLCHAR = 269,
    EOL = 270
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 632 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:355  */

 int intval;
 char * strval;
 double floatval;

#line 731 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE SAMlval;

int SAMparse (Extractor * state);

#endif /* !YY_SAM_HOME_VARTANIANMH_DEVEL_NCBI_VDB_LIBS_ALIGN_SAMEXTRACT_TOKENS_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 748 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   30

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  16
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  13
/* YYNRULES -- Number of rules.  */
#define YYNRULES  27
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  41

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   270

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   665,   665,   666,   670,   671,   676,   677,   678,   679,
     680,   681,   685,   691,   709,   722,   735,   746,   747,   750,
     758,   764,   770,   774,   778,   800,   801,   807
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "$undefined", "HEADER", "SEQUENCE",
  "READGROUP", "PROGRAM", "COMMENT", "TAG", "VALUE", "ALIGNVALUE", "QNAME",
  "COLON", "TAB", "CONTROLCHAR", "EOL", "$accept", "sam", "line",
  "comment", "header", "sequence", "program", "readgroup", "tagvaluelist",
  "tagvalue", "alignment", "avlist", "av", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270
};
# endif

#define YYPACT_NINF -12

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-12)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
     -12,     0,   -12,   -11,   -11,   -11,   -11,   -12,    -3,   -12,
     -12,   -12,   -12,    11,   -12,   -12,   -12,   -12,     4,   -11,
     -12,   -11,   -11,   -11,     6,    -3,   -12,   -12,     1,    -7,
     -12,   -12,   -12,   -12,   -12,    18,    16,   -12,   -12,    20,
     -12
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,    12,     0,     5,
       4,     3,     6,     0,     8,     9,    10,    11,     0,    13,
      17,    14,    16,    15,     0,    24,    25,     7,     0,     0,
      23,    18,    27,    26,    22,     0,     0,    21,    19,     0,
      20
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -12,   -12,   -12,   -12,   -12,   -12,   -12,   -12,    19,    -1,
     -12,   -12,     5
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    11,    12,    13,    14,    15,    16,    19,    20,
      17,    25,    26
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
       2,    36,    18,     3,     4,     5,     6,     7,    37,    34,
      24,     8,    28,    35,     9,    10,    32,    29,    31,    30,
      31,    31,    31,    21,    22,    23,    27,    38,    39,    40,
      33
};

static const yytype_uint8 yycheck[] =
{
       0,     8,    13,     3,     4,     5,     6,     7,    15,     8,
      13,    11,     8,    12,    14,    15,    10,    13,    19,    15,
      21,    22,    23,     4,     5,     6,    15,     9,    12,     9,
      25
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    17,     0,     3,     4,     5,     6,     7,    11,    14,
      15,    18,    19,    20,    21,    22,    23,    26,    13,    24,
      25,    24,    24,    24,    13,    27,    28,    15,     8,    13,
      15,    25,    10,    28,     8,    12,     8,    15,     9,    12,
       9
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    16,    17,    17,    18,    18,    18,    18,    18,    18,
      18,    18,    19,    20,    21,    22,    23,    24,    24,    25,
      25,    25,    25,    25,    26,    27,    27,    28
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     2,     1,     1,
       1,     1,     1,     2,     2,     2,     2,     1,     2,     4,
       5,     3,     3,     2,     2,     1,     2,     2
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (state, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, state); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, Extractor * state)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (state);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, Extractor * state)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, state);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, Extractor * state)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              , state);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, state); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, Extractor * state)
{
  YYUSE (yyvaluep);
  YYUSE (state);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (Extractor * state)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (state);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 5:
#line 671 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { ERR("CONTROLCHAR"); 
                   rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                   globstate->rc=rc;
                   return rc;
   }
#line 1851 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 6:
#line 676 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("comment"); }
#line 1857 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 7:
#line 677 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("header"); }
#line 1863 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 8:
#line 678 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("sequence"); }
#line 1869 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 9:
#line 679 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("program"); }
#line 1875 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 10:
#line 680 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("readgroup"); }
#line 1881 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 11:
#line 681 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("alignment"); }
#line 1887 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 12:
#line 685 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        mark_headers("CO");
    }
#line 1895 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 13:
#line 692 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("header tagvaluelist");
        check_required_tag(globstate->tags,"VN");
        if (!strcmp(globstate->tags,"SO ") &&
            !strcmp(globstate->tags,"GO "))
           WARN("Both SO and GO tags present");
        if (!(strcmp(globstate->tags,"SO ") ||
              strcmp(globstate->tags,"GO ")))
           WARN("neither SO or GO tags present");
        free(globstate->tags);
        globstate->tags=strdup("");

        mark_headers("HD");
    }
#line 1914 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 14:
#line 710 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("sequence");
        DBG(" sequences were: %s", globstate->seqnames);
        check_required_tag(globstate->tags,"SN");
        check_required_tag(globstate->tags,"LN");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("SQ");
    }
#line 1928 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 15:
#line 723 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("ids were: %s", globstate->ids);
        DBG("program");
        check_required_tag(globstate->tags,"ID");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("PG");
     }
#line 1941 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 16:
#line 736 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("readgroup");
        DBG("ids were: %s", globstate->ids);
        check_required_tag(globstate->tags,"ID");
        free(globstate->tags);
        globstate->tags=strdup("");
        mark_headers("RG");
     }
#line 1954 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 17:
#line 746 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG(" one tagvaluelist"); }
#line 1960 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 18:
#line 747 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG(" many tagvaluelist"); }
#line 1966 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 19:
#line 750 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("tagvalue:%s=%s", (yyvsp[-2].strval), (yyvsp[0].strval));
        const char * tag=(yyvsp[-2].strval);
        const char * value=(yyvsp[0].strval);
        process_tagvalue(tag,value);
        free((yyvsp[-2].strval));
        free((yyvsp[0].strval));
        }
#line 1979 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 20:
#line 758 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { 
        ERR("two tabs"); 
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        globstate->rc=rc;
        return rc;
  }
#line 1990 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 21:
#line 764 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { 
        ERR("empty tags"); 
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        globstate->rc=rc;
        return rc;
  }
#line 2001 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 22:
#line 770 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        const char * tag=(yyvsp[-1].strval);
        WARN("malformed TAG:VALUE 'TAB %s(NOT COLON)...'", tag);
        }
#line 2010 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 23:
#line 774 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { WARN("empty tags"); }
#line 2016 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 24:
#line 779 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG(" avlist qname:%s fields=%zu", (yyvsp[-1].strval), alignfields);
        alignfields=2;
        Alignment * align=myalloc(sizeof(Alignment));
        if (align==NULL) 
        {
            ERR("out of memory");
            rc_t rc=RC(rcAlign, rcRow,rcConstructing,rcMemory,rcExhausted);
            globstate->rc=rc;
            return rc;
        }
        align->read=globstate->read;
        align->cigar=globstate->cigar;
        align->rname=globstate->rname;
        align->pos=globstate->pos;
        VectorAppend(&globstate->alignments,NULL,align);
        free((yyvsp[-1].strval));
    }
#line 2039 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 25:
#line 800 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG(" one av"); }
#line 2045 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 26:
#line 801 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
           // TODO"bison: many avlist");
            }
#line 2053 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 27:
#line 808 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        const char * field=(yyvsp[0].strval);
        rc_t rc=process_align(field);
        globstate->rc=rc;
        free((yyvsp[0].strval));
    }
#line 2064 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;


#line 2068 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (state, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (state, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, state);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, state);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (state, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, state);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, state);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 816 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1906  */



 /* Epilogue */

