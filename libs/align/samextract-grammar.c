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

    #include <klib/rc.h>
    #include "samextract.h"
    #include <align/samextract-lib.h>
    int SAMlex(Extractor *);

#line 81 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:339  */

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
    VALUE = 263,
    QNAME = 264,
    FLAG = 265,
    RNAME = 266,
    POS = 267,
    MAPQ = 268,
    CIGAR = 269,
    RNEXT = 270,
    PNEXT = 271,
    TLEN = 272,
    SEQ = 273,
    QUAL = 274,
    OPTTAG = 275,
    OPTITAG = 276,
    OPTZTAG = 277,
    OPTBTAG = 278,
    OPTATYPE = 279,
    OPTITYPE = 280,
    OPTFTYPE = 281,
    OPTZTYPE = 282,
    OPTHTYPE = 283,
    OPTBTYPE = 284,
    OPTAVALUE = 285,
    OPTIVALUE = 286,
    OPTFVALUE = 287,
    OPTZVALUE = 288,
    OPTHVALUE = 289,
    OPTBVALUE = 290,
    HDVN = 291,
    HDSO = 292,
    HDGO = 293,
    RGID = 294,
    RGCN = 295,
    RGDS = 296,
    RGDT = 297,
    RGFO = 298,
    RGKS = 299,
    RGLB = 300,
    RGPG = 301,
    RGPI = 302,
    RGPL = 303,
    RGPM = 304,
    RGPU = 305,
    RGSM = 306,
    PLATFORM = 307,
    PGID = 308,
    PGPN = 309,
    PGCL = 310,
    PGPP = 311,
    PGDS = 312,
    PGVN = 313,
    SQSN = 314,
    SQLN = 315,
    SQAS = 316,
    SQM5 = 317,
    SQSP = 318,
    SQUR = 319,
    COLON = 320,
    TAB = 321,
    CONTROLCHAR = 322,
    EOL = 323
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 45 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:355  */

 char * strval;

#line 195 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE SAMlval;

int SAMparse (Extractor * state);

#endif /* !YY_SAM_HOME_VARTANIANMH_DEVEL_NCBI_VDB_LIBS_ALIGN_SAMEXTRACT_TOKENS_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 212 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:358  */

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
#define YYLAST   169

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  69
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  74
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  138

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   323

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
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   138,   138,   139,   143,   144,   148,   149,   150,   151,
     152,   153,   157,   161,   179,   180,   183,   187,   191,   195,
     199,   205,   223,   224,   228,   232,   242,   245,   250,   253,
     256,   260,   272,   273,   277,   281,   284,   287,   290,   293,
     296,   303,   316,   317,   320,   324,   327,   330,   333,   338,
     341,   344,   347,   350,   353,   357,   360,   363,   366,   371,
     375,   379,   385,   400,   417,   418,   422,   427,   432,   437,
     442,   447,   452,   457,   462
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "$undefined", "HEADER", "SEQUENCE",
  "READGROUP", "PROGRAM", "COMMENT", "VALUE", "QNAME", "FLAG", "RNAME",
  "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL",
  "OPTTAG", "OPTITAG", "OPTZTAG", "OPTBTAG", "OPTATYPE", "OPTITYPE",
  "OPTFTYPE", "OPTZTYPE", "OPTHTYPE", "OPTBTYPE", "OPTAVALUE", "OPTIVALUE",
  "OPTFVALUE", "OPTZVALUE", "OPTHVALUE", "OPTBVALUE", "HDVN", "HDSO",
  "HDGO", "RGID", "RGCN", "RGDS", "RGDT", "RGFO", "RGKS", "RGLB", "RGPG",
  "RGPI", "RGPL", "RGPM", "RGPU", "RGSM", "PLATFORM", "PGID", "PGPN",
  "PGCL", "PGPP", "PGDS", "PGVN", "SQSN", "SQLN", "SQAS", "SQM5", "SQSP",
  "SQUR", "COLON", "TAB", "CONTROLCHAR", "EOL", "$accept", "sam", "line",
  "comment", "header", "headerlist", "hdr", "sequence", "sequencelist",
  "sq", "program", "programlist", "pg", "readgroup", "readgrouplist", "rg",
  "alignment", "optlist", "opt", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323
};
# endif

#define YYPACT_NINF -48

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-48)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -48,     0,   -48,    53,    45,    37,    16,   -48,    -2,   -48,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,     6,    14,
      15,   -38,   -11,   -48,    21,    22,    23,    38,    46,    48,
     -48,    33,   -48,    52,    57,    58,    67,    90,    92,    94,
     102,   110,   112,    12,   113,   114,   115,   -47,    -7,   -48,
     -48,   116,   117,   118,   119,   120,   121,    -6,   -48,    42,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,    62,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48,
     -48,   122,   -48,   123,   124,   125,   126,   127,   128,   129,
     -10,    88,   106,   105,   104,   -48,    -5,   -48,   107,   108,
     103,   130,   109,   131,   133,   132,   134,   -48,   -48,   -48,
     -48,   -48,   -48,   -48,   -48,   -48,   -48,   -48
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,    12,     0,     5,
       4,     3,     6,     7,     8,     9,    10,    11,     0,     0,
       0,    20,     0,    14,     0,     0,     0,     0,     0,     0,
      30,     0,    22,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    42,
      40,     0,     0,     0,     0,     0,     0,     0,    32,     0,
      16,    17,    18,    19,    13,    15,    24,    25,    26,    27,
      28,    29,    21,    23,    58,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    54,    53,    55,    56,    57,    59,
      61,    41,    43,    34,    35,    36,    37,    38,    39,    31,
      33,     0,    60,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    62,     0,    64,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    63,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -48,   -48,   -48,   -48,   -48,   -48,   135,   -48,   -48,   136,
     -48,   -48,    84,   -48,   -48,    97,   -48,   -48,    31
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,    11,    12,    13,    22,    23,    14,    31,    32,
      15,    57,    58,    16,    48,    49,    17,   116,   117
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
       2,    33,    50,     3,     4,     5,     6,     7,    59,     8,
     111,   112,   113,   114,    60,   111,   112,   113,   114,    89,
      84,    90,    61,    62,    50,    18,    19,    20,    63,    66,
      67,    68,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    33,    69,    51,    52,    53,
      54,    55,    56,   101,    70,    21,    71,    64,   115,    47,
      74,    91,    99,   127,    85,    75,    76,     9,    10,    51,
      52,    53,    54,    55,    56,    77,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    18,
      19,    20,    24,    25,    26,    27,    28,    29,    78,    30,
      79,    72,    80,    47,    24,    25,    26,    27,    28,    29,
      81,    30,   118,   119,   120,   121,   122,   123,    82,    21,
      83,    86,    87,    88,    93,    94,    95,    96,    97,    98,
     102,   124,   125,   126,   103,   131,   104,   129,   105,   130,
     106,   100,   107,   133,   108,    92,   109,   128,   110,     0,
       0,     0,     0,     0,     0,     0,     0,    65,     0,     0,
       0,     0,     0,   132,   135,   136,   134,    73,     0,   137
};

static const yytype_int8 yycheck[] =
{
       0,     8,     8,     3,     4,     5,     6,     7,    10,     9,
      20,    21,    22,    23,     8,    20,    21,    22,    23,    66,
       8,    68,     8,     8,     8,    36,    37,    38,    66,     8,
       8,     8,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    50,    51,     8,     8,    53,    54,    55,
      56,    57,    58,    11,     8,    66,     8,    68,    68,    66,
       8,    68,    68,    68,    52,     8,     8,    67,    68,    53,
      54,    55,    56,    57,    58,     8,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,    36,
      37,    38,    59,    60,    61,    62,    63,    64,     8,    66,
       8,    68,     8,    66,    59,    60,    61,    62,    63,    64,
       8,    66,    24,    25,    26,    27,    28,    29,     8,    66,
       8,     8,     8,     8,     8,     8,     8,     8,     8,     8,
      68,    25,    27,    29,    12,    32,    13,    30,    14,    31,
      15,    57,    16,    34,    17,    48,    18,   116,    19,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    22,    -1,    -1,
      -1,    -1,    -1,    33,    31,    33,    35,    31,    -1,    35
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    70,     0,     3,     4,     5,     6,     7,     9,    67,
      68,    71,    72,    73,    76,    79,    82,    85,    36,    37,
      38,    66,    74,    75,    59,    60,    61,    62,    63,    64,
      66,    77,    78,     8,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    66,    83,    84,
       8,    53,    54,    55,    56,    57,    58,    80,    81,    10,
       8,     8,     8,    66,    68,    75,     8,     8,     8,     8,
       8,     8,    68,    78,     8,     8,     8,     8,     8,     8,
       8,     8,     8,     8,     8,    52,     8,     8,     8,    66,
      68,    68,    84,     8,     8,     8,     8,     8,     8,    68,
      81,    11,    68,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    68,    86,    87,    24,    25,
      26,    27,    28,    29,    25,    27,    29,    68,    87,    30,
      31,    32,    33,    34,    35,    31,    33,    35
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    69,    70,    70,    71,    71,    71,    71,    71,    71,
      71,    71,    72,    73,    74,    74,    75,    75,    75,    75,
      75,    76,    77,    77,    78,    78,    78,    78,    78,    78,
      78,    79,    80,    80,    81,    81,    81,    81,    81,    81,
      81,    82,    83,    83,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    84,    84,    84,    84,    84,
      84,    84,    85,    85,    86,    86,    87,    87,    87,    87,
      87,    87,    87,    87,    87
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     3,     1,     2,     2,     2,     2,     2,
       1,     3,     1,     2,     2,     2,     2,     2,     2,     2,
       1,     3,     1,     2,     2,     2,     2,     2,     2,     2,
       1,     3,     1,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       3,     2,    12,    13,     1,     2,     3,     3,     3,     3,
       3,     3,     3,     3,     3
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
#line 144 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { ERR("CONTROLCHAR");
                   rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
                   state->rc=rc;
   }
#line 1400 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 6:
#line 148 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("comment"); }
#line 1406 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 7:
#line 149 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("header done"); }
#line 1412 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 8:
#line 150 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("sequence"); }
#line 1418 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 9:
#line 151 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("program"); }
#line 1424 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 10:
#line 152 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("readgroup"); }
#line 1430 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 11:
#line 153 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("alignment"); }
#line 1436 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 12:
#line 157 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { mark_headers(state,"CO"); }
#line 1442 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 13:
#line 162 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
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
#line 1461 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 16:
#line 183 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->hashdvn=true;
        process_header(state,"HD","VN",(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1470 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 17:
#line 187 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->hashdso=true;
        process_header(state,"HD","SO",(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1479 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 18:
#line 191 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->hashdgo=true;
        process_header(state,"HD","GO",(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1488 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 19:
#line 195 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        ERR("two tabs"); /* TODO: Handle >2 tabs in a row */
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
#line 1497 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 20:
#line 199 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { WARN("empty tags"); }
#line 1503 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 21:
#line 206 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
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
#line 1523 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 24:
#line 228 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->hassqsn=true;
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1532 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 25:
#line 232 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        if (!inrange((yyvsp[0].strval),1,INT32_MAX))
        {
            ERR("SQ LN field not in range %s",(yyvsp[0].strval));
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        state->hassqln=true;
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1547 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 26:
#line 242 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1555 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 27:
#line 245 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        if (!ismd5((yyvsp[0].strval)))
            WARN("M5 value not followed by MD5");
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1565 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 28:
#line 250 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1573 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 29:
#line 253 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"SQ",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1581 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 30:
#line 256 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { WARN("Unexpected tab in sequence"); }
#line 1587 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 31:
#line 261 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("program");
        if (!state->haspgid)
        {
            ERR("ID tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }
        mark_headers(state,"PG"); }
#line 1601 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 34:
#line 277 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->haspgid=true;
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1610 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 35:
#line 281 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1618 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 36:
#line 284 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1626 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 37:
#line 287 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1634 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 38:
#line 290 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1642 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 39:
#line 293 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"PG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1650 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 40:
#line 296 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        WARN("Bogus value in PG:%s",(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1658 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 41:
#line 304 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("readgroup ");
        if (!state->hasrgid)
        {
            ERR("ID tag not seen in header");
            rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
            state->rc=rc;
        }

        mark_headers(state,"RG"); }
#line 1673 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 44:
#line 320 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        state->hasrgid=true;
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1682 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 45:
#line 324 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1690 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 46:
#line 327 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1698 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 47:
#line 330 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1706 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 48:
#line 333 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        if (!isfloworder((yyvsp[0].strval)))
            WARN("Flow order incorrec");
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1716 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 49:
#line 338 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1724 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 50:
#line 341 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1732 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 51:
#line 344 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1740 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 52:
#line 347 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1748 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 53:
#line 350 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1756 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 54:
#line 353 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        ERR("Invalid Platform %s", (yyvsp[0].strval));
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1765 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 55:
#line 357 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1773 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 56:
#line 360 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1781 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 57:
#line 363 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        process_header(state,"RG",(yyvsp[-1].strval),(yyvsp[0].strval));
        pool_free((yyvsp[0].strval)); }
#line 1789 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 58:
#line 366 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        WARN("Unknown readgroup (RG) tag:%s", (yyvsp[-1].strval));
        pool_free((yyvsp[-1].strval));
        pool_free((yyvsp[0].strval));
        }
#line 1799 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 59:
#line 371 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        ERR("two tabs"); /* TODO: Handle >2 tabs in a row */
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
#line 1808 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 60:
#line 375 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        ERR("empty tags");
        rc_t rc=RC(rcAlign,rcRow,rcParsing,rcData,rcInvalid);
        state->rc=rc; }
#line 1817 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 61:
#line 379 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { WARN("empty tags"); }
#line 1823 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 62:
#line 386 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("alignment record");
        process_alignment(state,(yyvsp[-11].strval),(yyvsp[-10].strval),(yyvsp[-9].strval),(yyvsp[-8].strval),(yyvsp[-7].strval),(yyvsp[-6].strval),(yyvsp[-5].strval),(yyvsp[-4].strval),(yyvsp[-3].strval),(yyvsp[-2].strval),(yyvsp[-1].strval));
        pool_free((yyvsp[-11].strval));
        pool_free((yyvsp[-10].strval));
        pool_free((yyvsp[-9].strval));
        pool_free((yyvsp[-8].strval));
        pool_free((yyvsp[-7].strval));
        pool_free((yyvsp[-6].strval));
        pool_free((yyvsp[-5].strval));
        pool_free((yyvsp[-4].strval));
        pool_free((yyvsp[-3].strval));
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[-1].strval)); }
#line 1842 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 63:
#line 401 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("alignment record with optional tags");
        process_alignment(state,(yyvsp[-12].strval),(yyvsp[-11].strval),(yyvsp[-10].strval),(yyvsp[-9].strval),(yyvsp[-8].strval),(yyvsp[-7].strval),(yyvsp[-6].strval),(yyvsp[-5].strval),(yyvsp[-4].strval),(yyvsp[-3].strval),(yyvsp[-2].strval));
        pool_free((yyvsp[-12].strval));
        pool_free((yyvsp[-11].strval));
        pool_free((yyvsp[-10].strval));
        pool_free((yyvsp[-9].strval));
        pool_free((yyvsp[-8].strval));
        pool_free((yyvsp[-7].strval));
        pool_free((yyvsp[-6].strval));
        pool_free((yyvsp[-5].strval));
        pool_free((yyvsp[-4].strval));
        pool_free((yyvsp[-3].strval));
        pool_free((yyvsp[-2].strval)); }
#line 1861 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 64:
#line 417 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG("opt"); }
#line 1867 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 65:
#line 418 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    { DBG(" opts"); }
#line 1873 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 66:
#line 423 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?AA");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1882 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 67:
#line 428 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?II");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1891 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 68:
#line 433 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?FF");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1900 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 69:
#line 438 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?ZZ");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1909 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 70:
#line 443 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?HH");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1918 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 71:
#line 448 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("?BB");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1927 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 72:
#line 453 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("III");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1936 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 73:
#line 458 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("ZZZ");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1945 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;

  case 74:
#line 463 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1646  */
    {
        DBG("BBB");
        pool_free((yyvsp[-2].strval));
        pool_free((yyvsp[0].strval)); }
#line 1954 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
    break;


#line 1958 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.c" /* yacc.c:1646  */
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
#line 469 "/home/vartanianmh/devel/ncbi-vdb/libs/align/samextract-grammar.y" /* yacc.c:1906  */


