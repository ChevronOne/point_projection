
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




#ifndef PPL_LERPER_HPP
#define PPL_LERPER_HPP


#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include "ppl_numeric_mth.hpp"
#include "boost_bar/progress.hpp"


namespace ppl { namespace LERPer

{


template<typename P_TYPE>
static const std::function<P_TYPE(const P_TYPE&)> basisF[ppl::cubic_points]={ 
            [](P_TYPE val) -> P_TYPE {return        std::pow(1.0-val, 3);},
            [](P_TYPE val) -> P_TYPE {return  3.0*std::pow(1.0-val, 2)*val;},
            [](P_TYPE val) -> P_TYPE {return        3.0*(1.0-val)*val*val;},
            [](P_TYPE val) -> P_TYPE {return              val*val*val;} };


template<typename P_TYPE> struct __state{
    uint64_t IND{0};
    P_TYPE ERR{0};
};

struct _sector{
    uint64_t f{0}, l{0}, len{0};
};


template<typename P_TYPE> PPL_FUNC_DECL void
_centripetal_par(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<P_TYPE>&,
        const ppl::LERPer::_sector&);


template<typename P_TYPE> PPL_FUNC_DECL ppl::vertex<P_TYPE> 
right_tan(const std::vector<ppl::vertex<P_TYPE>>&, const uint64_t&);


template<typename P_TYPE> PPL_FUNC_DECL ppl::vertex<P_TYPE> 
left_tan(const std::vector<ppl::vertex<P_TYPE>>&, const uint64_t&);



template<typename P_TYPE> PPL_FUNC_DECL ppl::LERPer::__state<P_TYPE> 
_attempt_to_fit(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LERPer::_sector&, 
        const ppl::LD&);


template<typename P_TYPE> ppl::LERPer::__state<P_TYPE>
fitSingle(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LD&);


template<typename P_TYPE> PPL_FUNC_DECL void 
_winnow(const std::vector<ppl::vertex<P_TYPE>>&,
        std::vector<uint64_t>&);


template<typename P_TYPE> PPL_FUNC_DECL ppl::LERPer::__state<P_TYPE> 
state_call(const std::vector<ppl::vertex<P_TYPE>>&, 
        const ppl::vertex<P_TYPE>* const,
        const ppl::LERPer::_sector&);


template<typename P_TYPE> ppl::LERPer::__state<P_TYPE>
extractB_path(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LD&);


template<typename P_TYPE> PPL_FUNC_DECL bool 
sector_tuneUP(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        ppl::LERPer::__state<P_TYPE>&,
        const ppl::LERPer::_sector&, 
        const ppl::LD&);

    
template<typename P_TYPE> PPL_FUNC_DECL void 
call_significant_fig_ascertain(const ppl::LD&);



#include "ppl_LERPer.inl"


} /* namespace LERPer */  



} /* namespace ppl */





#endif // PPL_LERPER_HPP


