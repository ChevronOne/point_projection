


/* 3D vector library */




#ifndef PPL_VERTEX_HPP
#define PPL_VERTEX_HPP


#include "ppl_headers.hpp"
#include <ostream>
#include <iomanip>
#include <sstream>
#include "ppl_debug.hpp"






namespace ppl
{

template<typename P_TYPE>
struct vertex{



    template<typename T>
    inline auto dot(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() * std::declval<T>());

    template<typename T>
    inline auto sqr_dist(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() + std::declval<T>());


    template<typename T>
    inline auto dist(const ppl::vertex<T>&) const 
        -> decltype(std::declval<P_TYPE>() + std::declval<T>());


    inline P_TYPE length(void) const;
    
    inline P_TYPE l1norm(void) const;

    inline P_TYPE l2norm(void) const; 

    inline std::string str(void) const;

    inline P_TYPE& operator[](const std::size_t);


    inline P_TYPE operator[](const std::size_t) const;


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator=(const T&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator+=(const ppl::vertex<T>&);

        
    template<typename T>
    inline ppl::vertex<P_TYPE>& operator-=(const ppl::vertex<T>&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator*=(const ppl::vertex<T>&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator/=(const ppl::vertex<T>&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator+=(const T&);

        
    template<typename T>
    inline ppl::vertex<P_TYPE>& operator-=(const T&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator*=(const T&);


    template<typename T>
    inline ppl::vertex<P_TYPE>& operator/=(const T&);



    vertex(const P_TYPE&, const P_TYPE&, const P_TYPE&);
    vertex(const ppl::vertex<P_TYPE>&);
    vertex(const P_TYPE&);

    vertex() = default;

    union{
		struct{ P_TYPE x, y, z; };
		struct{ P_TYPE r, g, b; };
		struct{ P_TYPE s, t, p; };
    };

}; // ppl::vertex




template<typename P_TYPE>
std::ostream& operator<<(std::ostream&, const ppl::vertex<P_TYPE>&);


template<typename P_TYPE>
inline std::string operator+(const std::string&, const ppl::vertex<P_TYPE>&);


template<typename P_TYPE>
inline std::string& operator+=
        (std::string& _str, 
        const ppl::vertex<P_TYPE>& vec);


template<typename T1, typename T2>
inline bool operator==(const ppl::vertex<T1>&, const ppl::vertex<T2>&);


template<typename T1, typename T2> 
inline bool operator==(const ppl::vertex<T1>&, const T2&);


template<typename T>
inline ppl::vertex<T> operator^(ppl::vertex<T>, std::size_t);


template<typename T>
inline ppl::vertex<T> operator-(const ppl::vertex<T>&);


template<typename T1, typename T2>
inline auto operator*(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;



template<typename T1, typename T2>
inline auto operator*(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator*(const T2&, const ppl::vertex<T1>&);



template<typename T1, typename T2>
inline auto operator+(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator+(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>;

template<typename T1, typename T2>
inline auto operator+(const T2&, const ppl::vertex<T1>&);



template<typename T1, typename T2>
inline auto operator-(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator-(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator-(const T2&, const ppl::vertex<T1>&) 
    -> ppl::vertex<decltype(std::declval<T2>() - std::declval<T1>())>;


template<typename T1, typename T2>
inline auto operator/(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator/(const ppl::vertex<T1>&, const T2&) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>;


template<typename T1, typename T2>
inline auto operator/(const T2&, const ppl::vertex<T1>&) 
    -> ppl::vertex<decltype(std::declval<T2>() / std::declval<T1>())>;


template<typename T>
inline ppl::vertex<T> sqrt(const ppl::vertex<T>&);


template<typename T>
inline ppl::vertex<T> normaize(const ppl::vertex<T>);


template<typename T1, typename T2>
inline auto cross(const ppl::vertex<T1>&, const ppl::vertex<T2>&) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>;

    
#include "ppl_vertex.inl"

} // namespace ppl





#endif  // PPL_VERTEX_HPP





