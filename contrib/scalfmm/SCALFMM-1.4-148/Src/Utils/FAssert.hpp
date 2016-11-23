// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FASSERT_HPP
#define FASSERT_HPP

#include <iostream>
#include <cassert>

#include "FGlobal.hpp"

// Clang warns about named variadic macros even though they are supported.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvariadic-macros"
#endif


/**
 * @brief The FError class
 * It is recommendede to use the macro:
 * FAssertLF( aTest , "some data ", "to ", plot);
 * Assertion are enabled or disabled during the compilation.
 * If disabled, the test instruction is still used (but the return will be optimized out by
 * the compiler).
 */
class FError {
protected:
    /**
     * @brief ErrPrint private method to end print
     */
    static void ErrPrint(){
        std::cerr << '\n';
    }

    /**
     * @brief ErrPrint private methdo to print
     */
    template<typename T, typename... Args>
    static void ErrPrint(const T& toPrint, Args... args){
        std::cerr << toPrint;
        ErrPrint( args... );
    }


public:
    //////////////////////////////////////////////////////////////
    // Should not be called directly
    //////////////////////////////////////////////////////////////

    /** Nothing to print */
    static void Print(){
    }

    /** One or more things to print */
    template<typename T, typename... Args>
    static void Print(const T& toPrint, Args... args){
        std::cerr << "[ERROR] ";
        ErrPrint( toPrint, args... );
    }
};

#ifdef SCALFMM_USE_ASSERT

//////////////////////////////////////////////////////////////
// Sp error activated
//////////////////////////////////////////////////////////////

#define FErrorAssertExit(TEST, args...) \
    if( !(TEST) ){ \
        FError::Print( args ); \
        throw std::exception(); \
    }


#else

//////////////////////////////////////////////////////////////
// Sp error desactivated
//////////////////////////////////////////////////////////////

#define FErrorAssertExit(TEST, args...) \
    if( !(TEST) ){}


#endif

//////////////////////////////////////////////////////////////
// Shortcut macro
//////////////////////////////////////////////////////////////

#define SPARSETD_ERROR_LINE " At line " , __LINE__ , "."
#define SPARSETD_ERROR_FILE " In file " , __FILE__ , "."

#define FAssert FErrorAssertExit

#define FAssertLF(args...) FAssert(args, SPARSETD_ERROR_LINE, SPARSETD_ERROR_FILE)

#endif //FASSERT_HPP


