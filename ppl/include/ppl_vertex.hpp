
/* 3D vector library */



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




#ifndef PPL_VERTEX_HPP
#define PPL_VERTEX_HPP



#include "ppl_debug.hpp"
#include <initializer_list>




namespace ppl
{

template<typename P_TYPE>
struct vertex{



    template<typename T>
    PPL_FUNC_DECL auto dot(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() * std::declval<T>());

    template<typename T>
    PPL_FUNC_DECL auto sqr_dist(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() + std::declval<T>());


    template<typename T>
    PPL_FUNC_DECL auto dist(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() + std::declval<T>());


    PPL_FUNC_DECL P_TYPE length(void) const;
    
    PPL_FUNC_DECL P_TYPE l1norm(void) const; 
    PPL_FUNC_DECL P_TYPE l2norm(void) const; 

    PPL_FUNC_DECL std::string str(void) const;

    PPL_FUNC_DECL P_TYPE& operator[](const std::size_t);


    PPL_FUNC_DECL P_TYPE operator[](const std::size_t) const;


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator=(const T&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator+=(const ppl::vertex<T>&);

        
    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator-=(const ppl::vertex<T>&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator*=(const ppl::vertex<T>&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator/=(const ppl::vertex<T>&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator+=(const T&);

        
    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator-=(const T&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator*=(const T&);


    template<typename T>
    PPL_FUNC_DECL ppl::vertex<P_TYPE>& operator/=(const T&);



    vertex(const P_TYPE&, const P_TYPE&, const P_TYPE&);
    vertex(const ppl::vertex<P_TYPE>&)= default;
    vertex(const P_TYPE&);
    vertex() = default;

    vertex(ppl::vertex<P_TYPE>&&)= default;
    vertex<P_TYPE>& operator=(const vertex<P_TYPE>&) = default;
    vertex<P_TYPE>& operator=(vertex<P_TYPE>&&) = default;


	union { P_TYPE x, r, s; };
	union { P_TYPE y, g, t; };
	union { P_TYPE z, b, p; };

}; // ppl::vertex



template<typename P_TYPE>
std::ostream& operator<<(std::ostream&, const ppl::vertex<P_TYPE>&);

template<typename P_TYPE>
PPL_FUNC_DECL std::string operator+(const std::string&, const ppl::vertex<P_TYPE>&);

template<typename P_TYPE>
PPL_FUNC_DECL std::string& operator+=
        (std::string& _str, 
        const ppl::vertex<P_TYPE>& vec);


template<typename T1, typename T2>
PPL_FUNC_DECL bool operator==(const ppl::vertex<T1>&, const ppl::vertex<T2>&);

template<typename T1, typename T2> 
PPL_FUNC_DECL bool operator==(const ppl::vertex<T1>&, const T2&);


template<typename T>
PPL_FUNC_DECL ppl::vertex<T> operator^(ppl::vertex<T>, std::size_t);

template<typename T>
PPL_FUNC_DECL ppl::vertex<T> operator-(const ppl::vertex<T>&);

template<typename T1, typename T2>
PPL_FUNC_DECL auto operator*(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;



template<typename T1, typename T2>
PPL_FUNC_DECL auto operator*(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator*(const T2&, const ppl::vertex<T1>&);



template<typename T1, typename T2>
PPL_FUNC_DECL auto operator+(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator+(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>;

template<typename T1, typename T2>
PPL_FUNC_DECL auto operator+(const T2&, const ppl::vertex<T1>&);



template<typename T1, typename T2>
PPL_FUNC_DECL auto operator-(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator-(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator-(const T2&, const ppl::vertex<T1>&) 
    -> ppl::vertex<decltype(std::declval<T2>() - std::declval<T1>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator/(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator/(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>;


template<typename T1, typename T2>
PPL_FUNC_DECL auto operator/(const T2&, const ppl::vertex<T1>&) 
    -> ppl::vertex<decltype(std::declval<T2>() / std::declval<T1>())>;


template<typename T>
PPL_FUNC_DECL ppl::vertex<T> sqrt(const ppl::vertex<T>&);

template<typename T>
PPL_FUNC_DECL ppl::vertex<T> normalize(const ppl::vertex<T>);

template<typename T1, typename T2>
PPL_FUNC_DECL auto cross(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;

    
#include "ppl_vertex.inl"

} // namespace ppl





#endif  // PPL_VERTEX_HPP





