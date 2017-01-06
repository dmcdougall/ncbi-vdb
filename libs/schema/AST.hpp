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

#ifndef _hpp_AST_
#define _hpp_AST_

#include "ParseTree.hpp"

#include <klib/symtab.h>

#include <vdb/schema.h>

struct KSymbol;

namespace ncbi
{
    namespace SchemaParser
    {
        class AST : public ParseTree
        {
        public:
            AST ( const Token* token );
            explicit AST (); // no token; a synthesized node

            virtual ~AST ();

            void AddNode ( AST * ); // allocate with new
            void AddNode ( const Token * );

            const AST* GetChild ( uint32_t idx ) const  { return static_cast < const AST * > ( ParseTree :: GetChild ( idx ) ); }
                  AST* GetChild ( uint32_t idx )        { return static_cast <       AST * > ( ParseTree :: GetChild ( idx ) ); }
        };

        class ASTBuilder
        {
        public:
            ASTBuilder ();
            ~ASTBuilder ();

            void DebugOn ();

            AST* Build ( const ParseTree& p_root );

            const VSchema * GetSchema () const { return m_schema; }

            void Declare ( const Token& p_ident, KSymbol* );
            void DeclareType ( const Token& p_ident, const KSymbol* p_super, const AST* p_dimension );

            KSymbol* Resolve ( const Token& p_ident );

            // error list is cleared by a call to Build
            uint32_t GetErrorCount() const { return VectorLength ( & m_errors ); }
            const char* GetErrorMessage ( uint32_t p_idx ) const { return ( const char * ) VectorGet ( & m_errors, p_idx ); }

            void ReportError ( const char* p_fmt, ... );

        private:
            bool Init();

        private:
            bool m_debug;

            VSchema*    m_intrinsic;
            VSchema*    m_schema;
            KSymTable   m_symtab;

            Vector      m_errors;
        };

        class AST_Schema : public AST
        {
        public:
            AST_Schema ( const Token*,
                         AST* decls, // NULL OK; is deleted here
                         unsigned int version = 1 );
            explicit AST_Schema ();

            virtual ~AST_Schema();

            void SetVersion ( const char* ); // version specified as n.m; checked for valid version number here

        private:
            uint32_t m_versMajor;
            uint32_t m_versMinor;
        };

        class AST_TypeDef : public AST
        {
        public:
            AST_TypeDef ( ASTBuilder &,
                          const Token*,
                          AST* baseType,
                          AST* newTypes );
            virtual ~AST_TypeDef();
        };

        class AST_ArrayDef : public AST
        {
        public:
            AST_ArrayDef ( const Token*,
                           AST* typeName,
                           AST* dimension );
            virtual ~AST_ArrayDef();
        };
    }
}

#endif
