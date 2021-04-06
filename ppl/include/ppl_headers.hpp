


//  Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
//  
//  This file is part of the Point Projection Library (ppl).
//  
//  Distributed under the terms of the GNU General Public License
//  as published by the Free Software Foundation; You should have
//  received a copy of the GNU General Public License.
//  If not, see <http://www.gnu.org/licenses/>.
//  
//  
//  This library is distributed in the hope that it will be useful, but WITHOUT
//  WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
//  WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
//  NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE
//  DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY,
//  WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
//  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. See the GNU
//  General Public License for more details.



/*
 * Copyright Abbas M.Murrey 2019-21
 *
 * Permission to use, copy, modify, distribute and sell this software
 * for any purpose is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both the copyright
 * notice and this permission notice appear in supporting documentation.  
 * I make no representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */



#ifndef PPL_HEADERS_HPP
#define PPL_HEADERS_HPP



#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string.h>
#include <limits>
#include <cmath>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <algorithm>
#include <functional>
#include <stdio.h>
#include <errno.h>

namespace ppl
{

        #define _USE_MATH_DEFINES
        #define CONST const
        #define PPL_COMPILER_DEPENDENT_INLINE inline
    
        typedef std::size_t UNS;
        typedef long double LD;
        typedef double S_PRECISION;


#ifdef _MSC_VER
        #define PPL_FORCEINLINE __forceinline   
#elif __GNUC__ | __clang__ | __MINGW32__ | __MINGW64__
        #define PPL_FORCEINLINE inline __attribute__((always_inline))
#else 
		#warning "'inlining' is not guaranteed!"
#define PPL_FORCEINLINE PPL_COMPILER_DEPENDENT_INLINE
#endif


#ifdef PPL_FORCE_INLINE 

#define PPL_FUNC_DECL PPL_FORCEINLINE

#else
#define PPL_FUNC_DECL PPL_COMPILER_DEPENDENT_INLINE 
#endif



CONST LD eps{ 1e-17L };
CONST UNS quintic{5};
CONST UNS quartic{4};
CONST UNS quintic_Coeffs{6};
CONST UNS cubic_points{4};
CONST UNS cubic{3};
CONST UNS cubic_Coeffs{4};
CONST UNS quadratic{2};
CONST UNS quadratic_Coeffs{3};
CONST UNS IT_OVR_FLW{3'000};
    
		
		
/*
 * User defined tolerance as an allowed margin of error for closest point.
 * Based on the precision of used data-type, the tolerance value can be
 * chosen as desired for any value in the open interval (0.0, 1.0) 'excluding 0.0 and 1.0'.
 * Regardless of the type of data being used, the tolerance has to be
 * always defined and initialized by a 'long double' value.
 */
    CONST LD TOLERANCE(0.00001L);
 



#define VERGENCE 0.001
#define NUM_OF(arr)  sizeof(arr) / sizeof(*arr)



} // namespace ppl



#ifdef PPL_EXTERNAL_TRACK_LOADING

#include <utility>
#include <fstream>

#if __has_include(<filesystem>)

#include <filesystem>
namespace fs = std::filesystem;



#elif __has_include(<experimental/filesystem>)

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem::v1;

#else
#error "the compiler dose not support C++ filesystem!"
#endif

//

#if __has_include(<sys/mman.h>)

// #include <experimental/string_view>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <iterator>

#else
#error "loading control points from external file is not supported on this system!"
#endif   


#endif



#ifdef PPL_CONCURRENCY

#if __has_include(<unistd.h>)
#include <unistd.h>
#endif

#include <stdarg.h>
#include <pthread.h>
#include <memory>

// #include <numeric>
#include "ppl_pthlib.hpp"
#include <cstdlib>


#include <tuple>

#if __has_include(<minwindef.h>)

#include <minwindef.h>

#endif

#if defined(_WIN32) || defined(WIN32)

#include<windows.h>
#include<Windows.h>

#endif

#endif



#endif  // PPL_HEADERS_HPP




