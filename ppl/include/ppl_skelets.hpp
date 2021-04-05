

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




#ifndef PPL_SKELETS_HPP
#define PPL_SKELETS_HPP



namespace ppl
{


template<typename P_TYPE> struct projection{
    ppl::vertex<P_TYPE> closest;
    uint64_t index;
    P_TYPE dist;
    P_TYPE parameter;
};

template<typename P_TYPE>
struct poly1d
{

    uint8_t d{ppl::quintic};
    P_TYPE coeffs[ppl::quintic_Coeffs];  
};

template<typename P_TYPE>
static const std::function<const ppl::vertex<P_TYPE>(const ppl::vertex<P_TYPE>* const)> 
    parametric_coeffs[ppl::cubic_points]={ 
        [](const ppl::vertex<P_TYPE>* const controlPs) -> const ppl::vertex<P_TYPE> {return   3*(controlPs[1]-controlPs[2]) + controlPs[3]-controlPs[0]; },
        [](const ppl::vertex<P_TYPE>* const controlPs) -> const ppl::vertex<P_TYPE> {return   3*(controlPs[0]+controlPs[2])-6*controlPs[1]; },
        [](const ppl::vertex<P_TYPE>* const controlPs) -> const ppl::vertex<P_TYPE> {return   3*(controlPs[1]-controlPs[0]); },
        [](const ppl::vertex<P_TYPE>* const controlPs) -> const ppl::vertex<P_TYPE> {return   *controlPs; } 
    };

template<typename P_TYPE> struct poly3d {
    std::size_t d{ppl::cubic};
    ppl::vertex<P_TYPE> coeffs[ppl::cubic_Coeffs]; 
};


template<typename P_TYPE>
struct deriv3d
{   std::size_t d{ppl::quadratic};
    ppl::vertex<P_TYPE> coeffs[ppl::quadratic_Coeffs]; 
};

template<typename P_TYPE> 
struct sturmSeq
{   sturmSeq()
    {   poly[0].coeffs[0]=1.0;
        poly[1].d=ppl::quartic;
        poly[1].coeffs[0]=0.0;
        poly[1].coeffs[1]=ppl::quintic;}
    std::size_t len{2};
    ppl::poly1d<P_TYPE> poly[ppl::quintic_Coeffs];
};

template<typename P_TYPE> 
struct leadingPolys
{
    leadingPolys(){
        poly[1].d = ppl::quartic;
    }
    ppl::poly1d<P_TYPE> poly[ppl::quintic_Coeffs];
};


template <typename P_TYPE> 
PPL_FUNC_DECL signed char __sign(P_TYPE val) {
    return (P_TYPE(0) < val) - (val < P_TYPE(0));
}


template<typename P_TYPE>
struct real_roots
{
    std::size_t num{0};
    P_TYPE zeros[ppl::quintic*20];

    PPL_FUNC_DECL void push(P_TYPE val){
        this->zeros[this->num] = val;
        this->num++;
    }

    PPL_FUNC_DECL void clear(void){
        this->num = 0;
    }
};


template<typename P_TYPE> PPL_FUNC_DECL P_TYPE 
exp_squaring(P_TYPE val, uint8_t __exp) // works only for exponents in interval [2,255]
{
    P_TYPE temp{1};
    while(__exp > 1){
        if(__exp%2 == 0) __exp/=2;
        else {
            temp*=val;
            __exp = (__exp-1)/2;
        }
        val*=val;
    }
    return val * temp;
}



PPL_FUNC_DECL uint64_t prec_call(CONST ppl::LD& val)
{
	if(val>=1)
		return 1;
	return (uint64_t)std::ceil(-std::log10(val)-ppl::eps);
}


template<class T> PPL_FUNC_DECL const char* __Tn(void) 
{
	if(strlen(typeid(T).name()) == 1)
	{
		switch(*typeid(T).name())
		{
			case 0x64: return "double";
			case 0x65: return "long double";
			case 0x66: return "float";
			default: return typeid(T).name();
		}

	}
	return typeid(T).name();
}




}  //  namespace ppl



#endif //  PPL_SKELETS_HPP











