


#ifndef PPL_LERPER_HPP
#define PPL_LERPER_HPP


#include "ppl_quintic_poly.hpp"
#include "ppl_verifications.hpp"
#include "ppl_vertex.hpp"
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <chrono>

#include "boost_bar/progress.hpp"


namespace ppl::LERPer
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


template<typename P_TYPE> inline void
_centripetal_par(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<P_TYPE>&,
        const ppl::LERPer::_sector&);


template<typename P_TYPE> inline ppl::vertex<P_TYPE> 
right_tan(const std::vector<ppl::vertex<P_TYPE>>&, const uint64_t&);


template<typename P_TYPE> inline ppl::vertex<P_TYPE> 
left_tan(const std::vector<ppl::vertex<P_TYPE>>&, const uint64_t&);



template<typename P_TYPE> inline ppl::LERPer::__state<P_TYPE> 
_attempt_to_fit(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LERPer::_sector&, 
        const ppl::LD&);


template<typename P_TYPE> inline ppl::LERPer::__state<P_TYPE>
fitSingle(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LD&);


template<typename P_TYPE> inline void 
_winnow(const std::vector<ppl::vertex<P_TYPE>>&,
        std::vector<uint64_t>&);


template<typename P_TYPE> inline ppl::LERPer::__state<P_TYPE> 
state_call(const std::vector<ppl::vertex<P_TYPE>>&, 
        const ppl::vertex<P_TYPE>* const,
        const ppl::LERPer::_sector&);


template<typename P_TYPE> inline ppl::LERPer::__state<P_TYPE>
extractB_path(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        const ppl::LD&);


template<typename P_TYPE> inline bool 
sector_tuneUP(const std::vector<ppl::vertex<P_TYPE>>&, 
        std::vector<ppl::vertex<P_TYPE>>&,
        ppl::LERPer::__state<P_TYPE>&,
        const ppl::LERPer::_sector&, 
        const ppl::LD&);

    
template<typename P_TYPE> inline void 
call_significant_fig_ascertain(const ppl::LD&);



#include "ppl_LERPer.inl"


} // namespace ppl::LERPer


#endif // PPL_LERPER_HPP


