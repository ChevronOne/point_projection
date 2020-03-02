
    
    
    
template<typename P_TYPE>
    template<typename T>
inline auto ppl::vertex<P_TYPE>::dot(const ppl::vertex<T>& vec) const
     -> decltype(std::declval<P_TYPE>() * std::declval<T>())
{
    return this->x*vec.x+this->y*vec.y+this->z*vec.z;
}

template<typename P_TYPE>
    template<typename T>
inline auto ppl::vertex<P_TYPE>::sqr_dist(const ppl::vertex<T>& vec) const 
    -> decltype(std::declval<P_TYPE>() + std::declval<T>())
{
    return (this->x-vec.x) * (this->x-vec.x)+ 
            (this->y-vec.y) * (this->y-vec.y)+
            (this->z-vec.z) * (this->z-vec.z);
}

template<typename P_TYPE>
    template<typename T>
inline auto ppl::vertex<P_TYPE>::dist(const ppl::vertex<T>& vec) const 
    -> decltype(std::declval<P_TYPE>() + std::declval<T>()) 
{
            
    return std::sqrt(std::pow(std::abs(this->x-vec.x),2)+std::pow(std::abs(this->y-vec.y),2)+std::pow(std::abs(this->z-vec.z),2));

}

template<typename P_TYPE>
inline P_TYPE ppl::vertex<P_TYPE>::length(void) const { 

    return std::sqrt(this->x*this->x + this->y*this->y + this->z*this->z);  
}
    
template<typename P_TYPE>
inline P_TYPE ppl::vertex<P_TYPE>::l1norm(void) const { 

    return std::abs(this->x) + std::abs(this->y) + std::abs(this->z);
}

template<typename P_TYPE>
inline P_TYPE ppl::vertex<P_TYPE>::l2norm(void) const {

    return this->length();
}


template<typename P_TYPE>
inline std::string ppl::vertex<P_TYPE>::str(void) const { 
    std::string str;
    return str+=*this;
}

template<typename P_TYPE>
inline P_TYPE& ppl::vertex<P_TYPE>::operator[](const std::size_t inx)
{
    if(inx == 0)
        return this->x;
    else if(inx==1)
        return this->y;
    else if(inx==2)
        return this->z;
    else
        _out_of_range("subscript out of range!");
}

template<typename P_TYPE>
inline P_TYPE ppl::vertex<P_TYPE>::operator[](const std::size_t inx) const
{
    if(inx == 0)
        return *this->x;
    else if(inx==1)
        return *this->y;
    else if(inx==2)
        return *this->z;
    else 
        _out_of_range("subscript out of range!");
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& ppl::vertex<P_TYPE>::operator=(const T& val) 
{
    this->x = this->y = this->z = val;
    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator+=(const ppl::vertex<T>& vec)
{
    this->x+=vec.x;
    this->y+=vec.y;
    this->z+=vec.z;
    return *this;
}
        
template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator-=(const ppl::vertex<T>& vec)
{
    this->x-=vec.x;
    this->y-=vec.y;
    this->z-=vec.z;

    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator*=(const ppl::vertex<T>& vec)
{
    this->x*=vec.x;
    this->y*=vec.y;
    this->z*=vec.z;

    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator/=(const ppl::vertex<T>& vec)
{
    this->x/=vec.x;
    this->y/=vec.y;
    this->z/=vec.z;

    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator+=(const T& val)
{
    this->x+=val;
    this->y+=val;
    this->z+=val;
    return *this;
}
        
template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator-=(const T& val)
{
    this->x-=val;
    this->y-=val;
    this->z-=val;

    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator*=(const T& val)
{
    this->x*=val;
    this->y*=val;
    this->z*=val;

    return *this;
}

template<typename P_TYPE>
    template<typename T>
inline ppl::vertex<P_TYPE>& 
    ppl::vertex<P_TYPE>::operator/=(const T& val)
{
    this->x/=val;
    this->y/=val;
    this->z/=val;

    return *this;
}

template<typename P_TYPE>
ppl::vertex<P_TYPE>::vertex(const P_TYPE& X, const P_TYPE& Y, const P_TYPE& Z)
    : x{X}, y{Y}, z{Z} {}

template<typename P_TYPE>
ppl::vertex<P_TYPE>::vertex(const ppl::vertex<P_TYPE>& vec)
    : x{vec.x}, y{vec.y},z{vec.z}{} 

template<typename P_TYPE>
ppl::vertex<P_TYPE>::vertex(const P_TYPE& val)
    : x{val}, y{val}, z{val} {} 





template<typename P_TYPE>
std::ostream& operator<<(std::ostream& out_s, const ppl::vertex<P_TYPE>& vec) 
{
    auto c_prec = out_s.precision();
    out_s.precision(std::numeric_limits<P_TYPE>::digits10);
    out_s << std::fixed 
        << vec.x << " "
        << vec.y << " "
        << vec.z;
    out_s.precision(c_prec);
    return out_s;
    
}

template<typename P_TYPE>
inline std::string operator+
    (const std::string& _str, const ppl::vertex<P_TYPE>& vec)
{
    std::stringstream str_stream;
    str_stream.precision(std::numeric_limits<P_TYPE>::digits10);
    str_stream << std::fixed <<
         " " << vec.x <<
         " " << vec.y << 
         " " << vec.z;

    return _str+str_stream.str();
    
}

template<typename P_TYPE>
inline std::string& operator+=
        (std::string& _str, 
        const ppl::vertex<P_TYPE>& vec){
    _str = _str+vec;
    return _str;
}


template<typename T1, typename T2>
inline bool operator==(const ppl::vertex<T1>& vec1, 
            const ppl::vertex<T2>& vec2)
{
    return { ((vec1.x==vec2.x) && (vec1.y==vec2.y) && (vec1.z==vec2.z))};
}

template<typename T1, typename T2> 
inline bool operator==(const ppl::vertex<T1>& vec, const T2& val)

{
    return { ((vec.x==val) && (vec.y==val) && (vec.z==val))};
}


template<typename T>
inline ppl::vertex<T> operator^(ppl::vertex<T> vec, std::size_t __exp)
{
    if(__exp == 2){
        return ppl::vertex<T>{vec.x*vec.x, vec.y*vec.y, vec.z*vec.z};
    }
    else
    {
        ppl::vertex<T> temp{1, 1, 1};
        while (__exp > 1){
            if (__exp % 2 == 0)
                __exp /= 2;
            else
            {
                temp *= vec;
                __exp = (__exp - 1) / 2;
            }
            vec *= vec;
        }
        return vec * temp;
    }
}

template<typename T>
inline ppl::vertex<T> operator-(const ppl::vertex<T>& vec)
{
    return {vec.x*(-1), vec.y*(-1), vec.z*(-1)};
}

template<typename T1, typename T2>
inline auto operator*(const ppl::vertex<T1>& vec1, const ppl::vertex<T2>& vec2) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>
{
    return {vec1.x*vec2.x, vec1.y*vec2.y, vec1.z*vec2.z};
}


template<typename T1, typename T2>
inline auto operator*(const ppl::vertex<T1>& vec, const T2& val) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>
{
    return {vec.x*val, vec.y*val, vec.z*val};
}

template<typename T1, typename T2>
inline auto operator*(const T2& val, const ppl::vertex<T1>& vec)
{
    return vec*val;
}


template<typename T1, typename T2>
inline auto operator+(const ppl::vertex<T1>& vec1, const ppl::vertex<T2>& vec2) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>
{
    return {vec1.x+vec2.x, vec1.y+vec2.y, vec1.z+vec2.z};
}

template<typename T1, typename T2>
inline auto operator+(const ppl::vertex<T1>& vec, const T2& val) 
    -> ppl::vertex<decltype(std::declval<T1>() + std::declval<T2>())>
{
    return {vec.x+val, vec.y+val, vec.z+val};
}

template<typename T1, typename T2>
inline auto operator+(const T2& val, const ppl::vertex<T1>& vec)
    {
        return vec+val;
    }


template<typename T1, typename T2>
inline auto operator-(const ppl::vertex<T1>& vec1, const ppl::vertex<T2>& vec2) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>
{
    return {vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z};
}

template<typename T1, typename T2>
inline auto operator-(const ppl::vertex<T1>& vec, const T2& val) 
    -> ppl::vertex<decltype(std::declval<T1>() - std::declval<T2>())>
{
    return {vec.x-val, vec.y-val, vec.z-val};
}

template<typename T1, typename T2>
inline auto operator-(const T2& val, const ppl::vertex<T1>& vec) 
    -> ppl::vertex<decltype(std::declval<T2>() - std::declval<T1>())>
{
    return {val-vec.x, val-vec.y, val-vec.z};
}

template<typename T1, typename T2>
inline auto operator/(const ppl::vertex<T1>& vec1, const ppl::vertex<T2>& vec2) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>
{
    return {vec1.x/vec2.x, vec1.y/vec2.y, vec1.z/vec2.z};
}

template<typename T1, typename T2>
inline auto operator/(const ppl::vertex<T1>& vec, const T2& val) 
    -> ppl::vertex<decltype(std::declval<T1>() / std::declval<T2>())>
{
    return {vec.x/val, vec.y/val, vec.z/val};
}

template<typename T1, typename T2>
inline auto operator/(const T2& val, const ppl::vertex<T1>& vec) 
    -> ppl::vertex<decltype(std::declval<T2>() / std::declval<T1>())>
{
    return {val/vec.x, val/vec.y, val/vec.z};
}

template<typename T>
inline ppl::vertex<T> sqrt(const ppl::vertex<T>& vec){
    return { std::sqrt(vec.x), std::sqrt(vec.y), std::sqrt(vec.z)};
}

template<typename T>
inline ppl::vertex<T> normaize(const ppl::vertex<T> vec){
    return vec/std::sqrt(vec.dot(vec));
}

template<typename T1, typename T2>
inline auto cross(const ppl::vertex<T1>& vec1, const ppl::vertex<T2>& vec2) 
    -> ppl::vertex<decltype(std::declval<T1>() * std::declval<T2>())>
{
    return {vec1.y * vec2.z - vec2.y * vec1.z,
			vec1.z * vec2.x - vec2.z * vec1.x,
			vec1.x * vec2.y - vec2.x * vec1.y};
}





