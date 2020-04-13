

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
 * Copyright Abbas M.Murrey 2019-20
 *
 * Permission to use, copy, modify, distribute and sell this software
 * for any purpose is hereby granted without fee, provided that the
 * above copyright notice appear in all copies and that both the copyright
 * notice and this permission notice appear in supporting documentation.  
 * I make no representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */



#ifndef PPL_NUMERIC_MTH_HPP
#define PPL_NUMERIC_MTH_HPP

#include "ppl_skelets.hpp"

#include <initializer_list>
#include <algorithm>
#include <math.h>

namespace ppl{



template<typename P_TYPE>
class cubic_path
{    

    ppl::vertex<P_TYPE> (*splines)[ppl::cubic_points]{nullptr}; 

    ppl::vertex<P_TYPE>* control_points{nullptr};

    ppl::poly1d<P_TYPE>* polys{nullptr};

    ppl::poly3d<P_TYPE>* parametric{nullptr};
    ppl::deriv3d<P_TYPE>* deriv{nullptr};
    
    uint64_t poly_num{0};
    uint64_t points_num{0};


    const P_TYPE TOLERZ{
        static_cast<P_TYPE>(ppl::TOLERANCE<ppl::eps? ppl::eps:ppl::TOLERANCE)
    };

    const std::size_t NEWTON_THRES{
        (std::size_t)std::ceil( 
            std::log10(std::ceil(
                -std::log10(TOLERZ)) ) / std::log10(2.0)) + 1 };

    const P_TYPE upperB{static_cast<P_TYPE>(1.0 - TOLERZ)};
    const P_TYPE lowerB{TOLERZ};

    const uint32_t min_depth{ppl::quintic};


    const std::function<const P_TYPE(const uint64_t&)> 
        object_poly_coeffs[ppl::quintic_Coeffs]={ 
            [this](const uint64_t& i) -> const P_TYPE {return  -3 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[0]); } ,
            [this](const uint64_t& i) -> const P_TYPE {return (-5 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[1]))  / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return (-4 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[2])
                                                               -2 *  parametric[i].coeffs[1].dot(parametric[i].coeffs[1]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return  -3 * (parametric[i].coeffs[0].dot(parametric[i].coeffs[3])
                                                                   + parametric[i].coeffs[1].dot(parametric[i].coeffs[2]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return (-2 *  parametric[i].coeffs[1].dot(parametric[i].coeffs[3])
                                                                   - parametric[i].coeffs[2].dot(parametric[i].coeffs[2]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return       -parametric[i].coeffs[2].dot(parametric[i].coeffs[3])   / polys[i].coeffs[0]; }
        };


    void extract_poly(const uint64_t& _stride, const ppl::vertex<P_TYPE>* const points)
    {

        uint32_t i{0};
        for(i=0; i<ppl::cubic_points; ++i)
            parametric[_stride].coeffs[i] = ppl::parametric_coeffs<P_TYPE>[i](points);

        if(parametric[_stride].coeffs[0] == 0.0) 
            throw_arg_exception(_stride, points);

        for(i=0; i<ppl::cubic; ++i)
            deriv[_stride].coeffs[i] = (ppl::cubic-i) * parametric[_stride].coeffs[i];

        for(i=0;i<ppl::quintic_Coeffs; ++i)
            polys[_stride].coeffs[i] = object_poly_coeffs[i](_stride);

    }

    void throw_arg_exception(const std::size_t ind, 
            const ppl::vertex<P_TYPE>* const points) const 
    {

        // out << "9-0-0870";
        std::string err = "control points are on top of each other!\nstarting from the "
                    + std::to_string((ind * ppl::cubic)+1) 
                    + "'th control point:\n";

        for(std::size_t i{0}; i < ppl::cubic_points; ++i ){
            err+=std::to_string(i + (ind * ppl::cubic) +1);
            err += *(points+i);
        }

        std::cerr<< err << "\n";
        ppl_invalid_argument("invalid argument");
    }


    void call_significant_fig_ascertain(void) {

/*	To avoid 'constant expression' warning generated by MSVC! 
 *	This can be replaced by a single 'if constexpr', but since
 *	not all distributions support C++17 as default, this could
 *	be a temporary solution until the compilers catch up.
 */
#ifdef _MSC_VER
__pragma(warning(push))
__pragma(warning(disable:4127))
        if (!(ppl::TOLERANCE>0) || ppl::TOLERANCE>=1){
__pragma(warning(pop))
#else
		if (!(ppl::TOLERANCE > 0) || ppl::TOLERANCE >= 1){
#endif
			ppl_out_of_range("invalid tolerance value! the value of tolerance has to be in the open interval (0.0, 1.0) excluding 0.0 and 1.0 \n");
        }

		std::streamsize t_prec = static_cast<std::streamsize>(-std::log10(ppl::TOLERANCE));
        if( ppl::prec_call(ppl::TOLERANCE) > std::numeric_limits< P_TYPE>::digits10 ){ 
            
			std::cerr.precision(t_prec + 1);
            std::cerr << "the tolerance value of "
                      << std::fixed << ppl::TOLERANCE 
                      << " is out of range for precision of type <data type " 
                      << ppl::__Tn<P_TYPE>() << "> \n";

            ppl_logic_error("incompatible arguments");
        }
        
    }

    PPL_FUNC_DECL bool newton_mth(const ppl::leadingPolys<P_TYPE>& lead_polys,
                        P_TYPE val, 
                        const P_TYPE& a, const P_TYPE& b, 
                        ppl::real_roots<P_TYPE>& roots) const  
    {

        P_TYPE polyEvalu, derivEvalu;

        for (std::size_t i{0};;++i)
        {   
            polyEvalu = poly1d_solve_for(lead_polys.poly[0], val); 

            if ( std::abs( polyEvalu ) <= TOLERZ){
                roots.push(val);
                return 0;
            }

            derivEvalu = poly1d_solve_for(lead_polys.poly[1], val);

            if ( derivEvalu == 0.0 || i > NEWTON_THRES)  //  <<<<<<<<<<<<<<<< NEWTON'S METHOD FAILED!!  
                return 1;  // >>>>>> throw local maximum/minimum || iteration overflow!

            val = val - ( polyEvalu / derivEvalu );

            if(val < a || val > b) //     <<<<<<<<<<<<<  NEWTON'S METHOD FAILED!!
                return 1;  //  >>>>> throw wrong root
        }

    }

    template<typename sturm_t>
    PPL_FUNC_DECL uint8_t _signs_at
                (const ppl::sturmSeq<sturm_t>& seq,
                const P_TYPE &val_t) const {
        
        sturm_t val = static_cast<sturm_t>(val_t);
		uint8_t alters{0};
		bool priv, curr;
        priv = std::signbit(poly1d_solve_for(seq.poly[0],val));
    
        for (uint8_t i{1}; i < seq.len; ++i){
            curr =  std::signbit(poly1d_solve_for(seq.poly[i],val));
            alters+=(priv^curr);
            priv = curr;
        }

        return alters;
    }

    
    template<typename sturm_t>
    PPL_FUNC_DECL uint8_t roots_in_interval
                (const ppl::sturmSeq<sturm_t>& seq,
                const P_TYPE& a, const P_TYPE& b) const {
        return (uint8_t)std::abs(_signs_at(seq, a)-_signs_at(seq, b));
    }

    template<typename sturm_t>
    inline uint8_t __split(const ppl::sturmSeq<sturm_t>& seq, 
                P_TYPE a, P_TYPE b, 
                const uint8_t& _rN, 
                ppl::real_roots<P_TYPE>& roots,
                const ppl::leadingPolys<P_TYPE>& lead_polys,
                uint32_t curr_depth) const
    {
        if ((b - a) <= TOLERZ){
            roots.push( (a + b) / 2 ); 
            return _rN;
        }

        P_TYPE m_value = (a+b) /2;

        if (_rN == 1)
        {


            P_TYPE r_evalu{poly1d_solve_for(lead_polys.poly[0], b)};
            
            if (poly1d_solve_for(lead_polys.poly[0], a) < 0.0 && r_evalu > 0.0)
            {
                /*
                    Newton's method is extremely fast to find a root, 
                    but if it FOR VERY RARE SITUATION failed to find a root in a certain number of iterations, 
                    then more likely it'll not find a root at all, or it could tend toward a wrong root!
                    For some situations such as oscillating sequence it's a must to change 
                    the initial value by shrinking the interval using Bisection method.
                    
                */

               for(;curr_depth<min_depth;++curr_depth){
                   if (ppl::__sign(poly1d_solve_for(lead_polys.poly[0], 
                            m_value)) == ppl::__sign(r_evalu)){

                        b = m_value;
                        r_evalu = poly1d_solve_for(lead_polys.poly[0], b);
                    }
                    else a = m_value;

                    m_value = (a+b) / 2;
               }

                for(;newton_mth(lead_polys, m_value, a, b, roots);){

                    if (ppl::__sign(poly1d_solve_for(lead_polys.poly[0], 
                                    m_value)) == ppl::__sign(r_evalu)){
                        b = m_value;
                        r_evalu = poly1d_solve_for(lead_polys.poly[0], b);
                    }
                    else a = m_value;

                    m_value = (a+b) /2;

                    if ((b - a) <= TOLERZ){
                        roots.push( m_value ); 
                        return 1;
                    }
                }
                return 1;

            } else return 1;


        }else {

            uint8_t rootsN = roots_in_interval(seq, m_value + TOLERZ, b);
            ++curr_depth;
            if( rootsN != 0 ){
                rootsN = __split(seq, m_value + TOLERZ,
                            b, rootsN, roots, lead_polys, curr_depth);
            }

            if(rootsN != _rN) 
                return rootsN + __split(seq, a, 
                        m_value, _rN -rootsN, roots, lead_polys, curr_depth);
            else return rootsN;
        }
    }


    template<typename sturm_t>
    PPL_FUNC_DECL bool poly_long_division
            (const ppl::poly1d<sturm_t>& numer,
            const ppl::poly1d<sturm_t>& denom, 
            ppl::poly1d<sturm_t>& rem) const
    {

        memcpy(rem.coeffs, numer.coeffs, 
                sizeof(sturm_t) * ppl::quintic_Coeffs);
        rem.d = numer.d;

        for ( ;denom.d<=rem.d; )
        {
			sturm_t QUOT{ rem.coeffs[ppl::quintic - rem.d]
					/ denom.coeffs[ppl::quintic - denom.d] };

            rem.coeffs[ppl::quintic - rem.d] = 0.0;
			std::size_t i{ ppl::quintic - rem.d + 1 };
			std::size_t j{ i + denom.d };
            for (;i < j; ++i)
                rem.coeffs[i] = rem.coeffs[i]
                    - (QUOT * denom.coeffs[rem.d - denom.d + i ] );

            --rem.d;

            for (i = ppl::quintic - rem.d ; i < ppl::quintic_Coeffs; ++i)
                if (rem.coeffs[i] == 0.0){
                    if(rem.d == 0)
                        return 0;
                    else --rem.d;

                }else break;
        }

        std::for_each(rem.coeffs+(ppl::quintic-rem.d),
                    rem.coeffs+ppl::quintic_Coeffs, 
                    [](sturm_t& coeff){ coeff = -coeff;});

        return 1;
    }

    template<typename sturm_t>
    PPL_FUNC_DECL void construct_sturmPolys(ppl::sturmSeq<sturm_t>& seq) const {
        
        std::size_t i;
        for (i = 1; i < ppl::quintic; ++i) // first derivative
            seq.poly[1].coeffs[i + 1] = (seq.poly[0].d - i) * seq.poly[0].coeffs[i];

        for (i = 1; seq.poly[i].d > 0; ++i)
            if(poly_long_division(seq.poly[i-1], seq.poly[i], seq.poly[i+1]))
                seq.len++;
    }

    template<typename sturm_t>
    PPL_FUNC_DECL sturm_t poly1d_solve_for
            (const ppl::poly1d<sturm_t>& poly,
            const sturm_t& val) const
    {

        sturm_t result = poly.coeffs[ppl::quintic-poly.d];

        std::for_each(poly.coeffs+(ppl::quintic-poly.d + 1), 
                        poly.coeffs+ppl::quintic_Coeffs,
                            [&](const sturm_t& coeff){ result = result*val + coeff;}
                        );

        return result;
    }

    PPL_FUNC_DECL ppl::vertex<P_TYPE> poly3d_solve_for
            (const ppl::poly3d<P_TYPE>& poly, 
            const P_TYPE& val) const{

        return val*(val*(val*poly.coeffs[0]
                +poly.coeffs[1])+poly.coeffs[2]) + poly.coeffs[3];
    }

    template<typename alters_t>
    PPL_FUNC_DECL uint8_t __alters_in_coeffs
                (alters_t const * const _poly) const {
        uint16_t alters{0};
        bool priv, curr;
        priv = std::signbit(_poly[0]);
    
        for (uint8_t i{1}; i < ppl::quintic_Coeffs; ++i){
            curr =  std::signbit(_poly[i]);
            alters+=(priv^curr);
            priv = curr;
        }

        return alters;
    }


    PPL_FUNC_DECL void _call_projection(ppl::vertex<P_TYPE> const * const p, 
                ppl::projection<P_TYPE> * const point_projection) const
    {

        ppl_assert__(poly_num>0, 
            "closest point was called on empty data! did you forget to load your data?\n");

        point_projection->closest = splines[0][0];
		point_projection->index = 0;
		point_projection->parameter = 0;
        P_TYPE min_dist{(*p).sqr_dist( splines[0][0] )}, curr_dist;
    
        for (std::size_t i{0}; i < poly_num; ++i)
        {

            curr_dist = (*p).sqr_dist(splines[i][ppl::cubic]);
            if (min_dist > curr_dist){
                point_projection->closest = splines[i][ppl::cubic];
                point_projection->index = i;
                point_projection->parameter = static_cast<P_TYPE>(1);
                min_dist = curr_dist;
            }

            ppl::sturmSeq<ppl::S_PRECISION> seq;
            for(uint8_t j{1}; j<ppl::quintic_Coeffs; ++j)
                if(j<ppl::cubic)
                    seq.poly[0].coeffs[j] = polys[i].coeffs[j];
                else
                    seq.poly[0].coeffs[j] = polys[i].coeffs[j] 
                        + (deriv[i].coeffs[j-ppl::cubic].dot(*p) 
                        / polys[i].coeffs[0]);

            construct_sturmPolys(seq);

            uint8_t _rN{roots_in_interval(seq, lowerB, upperB)};

            if (_rN != 0){
                ppl::real_roots<P_TYPE> roots;
                ppl::leadingPolys<P_TYPE> lead_polys;
                for(std::size_t f{0};f<ppl::quintic_Coeffs; ++f)
                    lead_polys.poly[0].coeffs[f] = static_cast<P_TYPE>(seq.poly[0].coeffs[f]);
                for(std::size_t f{1};f<ppl::quintic_Coeffs; ++f)
                    lead_polys.poly[1].coeffs[f] = static_cast<P_TYPE>(seq.poly[1].coeffs[f]);
            
                __split(seq, lowerB, upperB, _rN, roots, lead_polys, 1);

                for (std::size_t j{0}; j < roots.num; ++j){
                    curr_dist = (*p).sqr_dist(poly3d_solve_for(parametric[i], 
                                                roots.zeros[j]));
                    if (min_dist > curr_dist){
                        point_projection->index = i;
                        point_projection->parameter = roots.zeros[j];
						min_dist = curr_dist;
                    }
                }
            }  

        }
    }


    void cleanUp(void)
    {
        __freem( splines);
        __freem( control_points);
        __freem( polys);
        __freem( parametric);
        __freem( deriv);

        poly_num = points_num = 0;

    }

    template< typename T > 
    PPL_FUNC_DECL void __freem(T* &_alloc){
        if(_alloc != nullptr){
            delete[] _alloc;
            _alloc = nullptr;
        }
    }

public:    
    // cubic_path() = default;
    cubic_path(): 
            splines{nullptr}, control_points{nullptr},   
            polys{nullptr}, parametric{nullptr}, 
            deriv{nullptr},
            poly_num{0}, points_num{0}, min_depth{ppl::quintic}
    {
        ppl_assert__(std::numeric_limits<P_TYPE>::is_iec559, 
            "instantiation of ppl::cubic_path can only be with floating-point types!\n");
    }

    cubic_path(const ppl::vertex<P_TYPE>* const points,
            const uint64_t& _size): min_depth{ppl::quintic} {
        ppl_assert__(std::numeric_limits<P_TYPE>::is_iec559, 
            "instantiation of ppl::cubic_path can only be with floating-point types!\n");
        routing(points, _size);
    }
       
    virtual ~cubic_path() { cleanUp(); }


#ifdef __CONCURRENCY__

    void closest_point(ppl::vertex<P_TYPE> const * const p, 
            ppl::projection<P_TYPE> * const projection_ptr) const
    {
        _call_projection(p, projection_ptr);
        if(projection_ptr->parameter != static_cast<P_TYPE>(1) && projection_ptr->parameter != 0){
            projection_ptr->closest = poly3d_solve_for(parametric[projection_ptr->index], 
                                                    projection_ptr->parameter);
        }

        projection_ptr->dist = projection_ptr->closest.dist(*p);
    }

#endif

    void routing(const ppl::vertex<P_TYPE>* const points, 
                    const uint64_t& _size)
    {      
        ppl_assert__( (_size -1)%ppl::cubic == 0 && _size > ppl::cubic, 
                    "incompatible number of control points!");
        call_significant_fig_ascertain();

        if(points_num != 0)  cleanUp();

        points_num = _size;
        poly_num = (_size-1)/ppl::cubic; 

        
        control_points = new ppl::vertex<P_TYPE>[_size];
        memcpy(control_points, points, sizeof(ppl::vertex<P_TYPE>)*_size); 

        splines = new ppl::vertex<P_TYPE>[poly_num] [ppl::cubic_points];
        
        std::size_t term = points_num-ppl::cubic;
        uint64_t i{0}, j{0};
        for(i=0, j=0; i < term; ++j, i+=ppl::cubic){
            memcpy(splines[j], points+i, 
                sizeof(ppl::vertex<P_TYPE>)*ppl::cubic_points);
        }

        polys = new ppl::poly1d<P_TYPE>[poly_num];

        parametric = new ppl::poly3d<P_TYPE>[poly_num];
        deriv = new ppl::deriv3d<P_TYPE>[poly_num];

        for(i=0, j=0; i < term; ++j, i+=ppl::cubic)
            extract_poly(j, points+i);

    }

    ppl::projection<P_TYPE> 
    closest_point(ppl::vertex<P_TYPE> const * const p) const

    {
        ppl::projection<P_TYPE> point_projection;
        _call_projection(p, &point_projection);
        if(point_projection.parameter != static_cast<P_TYPE>(1) && point_projection.parameter != 0){
            point_projection.closest = poly3d_solve_for(parametric[point_projection.index], 
                                        point_projection.parameter);
        }

        point_projection.dist = point_projection.closest.dist(*p);
        return point_projection;
    
    }

};
    
} // namespace ppl

#endif //  PPL_NUMERIC_MTH_HPP


