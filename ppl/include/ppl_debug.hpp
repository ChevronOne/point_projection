


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




#ifndef PPL_DEBUG_HPP
#define PPL_DEBUG_HPP

#include "ppl_headers.hpp"

namespace ppl
{


#if defined(_MSC_VER)

#define ppl_BUF_SIZE 0x01F4

char ppl_MSG_BUF[ppl_BUF_SIZE];

#define ppl_error_msg [&_BUF=ppl_MSG_BUF]() -> const char* { strerror_s(_BUF, ppl_BUF_SIZE, errno); return _BUF; }

#define ppl_reset_errNo() (errno == 0 ? "None" : ppl_error_msg())

#else

#define ppl_reset_errNo() (errno == 0 ? "None" : strerror(errno))

#endif

#define ppl_log_err(M, ...) fprintf(stderr,\
        "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__,\
        ppl_reset_errNo(), ##__VA_ARGS__)

#define ppl_assert__(A, M, ...) if(!(A)) {\
    ppl_log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define ppl_out_of_range(M, ...) {\
    ppl_log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define ppl_logic_error(M, ...) {\
    ppl_log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define ppl_invalid_argument(M, ...) {\
    ppl_log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define ppl_runtime_error(M, ...) {\
    ppl_log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _incompatible_platform(M, ...) {\
    std::cerr << "[ERROR] " << "errno:"<< ppl_reset_errNo()\
              << ",  "<< M << "\n"; errno=0; exit(EXIT_FAILURE); }


} // namespace ppl

#endif // PPL_DEBUG_HPP




