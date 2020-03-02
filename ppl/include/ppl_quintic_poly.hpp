


#ifndef QUINTIC_POLY_HPP
#define QUINTIC_POLY_HPP


#if defined(__CONCURRENCY__) &&  !defined(_WIN32) && !defined(WIN32) && !defined(__linux__)
#error "multithreading support is not supported on this system!"

#endif

#if defined(__EXTERNAL_TRACK_LOADING__) && !defined(__linux__)

#error "loading control points from external file is not supported on this system!"
#endif


#include "ppl_headers.hpp"
#include "ppl_vertex.hpp"
#include <initializer_list>
#include <algorithm>
#include <math.h>
#include "ppl_debug.hpp"
#include "ppl_verifications.hpp"





namespace ppl{
    

template<typename P_TYPE> struct projection{
    ppl::vertex<P_TYPE> closest;
    uint64_t index;
    P_TYPE dist;
    P_TYPE parameter;
};

template<typename P_TYPE>
struct poly1d
{

    std::size_t d{ppl::quintic};
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

template<typename P_TYPE> struct sturmSeq
{   sturmSeq()
    {   poly[0].coeffs[0]=1.0;
        poly[1].d=ppl::quartic;
        poly[1].coeffs[0]=0.0;
        poly[1].coeffs[1]=ppl::quintic;}
    std::size_t len{2};
    ppl::poly1d<P_TYPE> poly[ppl::quintic_Coeffs];
};


template <typename P_TYPE> inline signed char __sign(P_TYPE val) {
    return (P_TYPE(0) < val) - (val < P_TYPE(0));
}


template<typename P_TYPE>
struct real_roots
{
    std::size_t num{0};
    P_TYPE zeros[ppl::quintic*20];

    inline void push(P_TYPE val){
        this->zeros[this->num] = val;
        this->num++;
    }

    inline void clear(void){
        this->num = 0;
    }
};


template<typename P_TYPE> inline P_TYPE exp_squaring(P_TYPE val, uint8_t __exp) // works only for exponents in interval [2,255]
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


template<typename P_TYPE>
class cubic_path
{    

    ppl::vertex<P_TYPE> (*splines)[cubic_points]{nullptr}; 

    ppl::vertex<P_TYPE>* control_points{nullptr};

    ppl::poly1d<P_TYPE>* polys{nullptr};

    ppl::poly3d<P_TYPE>* parametric{nullptr};
    ppl::deriv3d<P_TYPE>* deriv{nullptr};
    
    uint64_t poly_num{0};
    uint64_t points_num{0};

    
    ppl::projection<P_TYPE> point_projection;
    P_TYPE min_dist, curr_dist;

    const P_TYPE TOLERZ{ppl::TOLERANCE<ppl::eps? ppl::eps:ppl::TOLERANCE};


    const P_TYPE upperB{static_cast<P_TYPE>(1.0 - TOLERZ)};
    const P_TYPE lowerB{static_cast<P_TYPE>(TOLERZ)};

    const std::size_t NEWTON_THRES{(std::size_t)std::ceil( std::log10(std::ceil(-std::log10(TOLERZ)) ) /log10(2) )+1};



    const std::function<const P_TYPE(const uint64_t&)> 
        object_poly_coeffs[ppl::quintic_Coeffs]={ 
            [this](const uint64_t& i) -> const P_TYPE {return  -3 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[0]); } ,
            [this](const uint64_t& i) -> const P_TYPE {return (-5 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[1])) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return (-4 *  parametric[i].coeffs[0].dot(parametric[i].coeffs[2])- 2*parametric[i].coeffs[1].dot(parametric[i].coeffs[1]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return  -3 * (parametric[i].coeffs[0].dot(parametric[i].coeffs[3])+   parametric[i].coeffs[1].dot(parametric[i].coeffs[2]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return (-2 *  parametric[i].coeffs[1].dot(parametric[i].coeffs[3])-   parametric[i].coeffs[2].dot(parametric[i].coeffs[2]) ) / polys[i].coeffs[0]; },
            [this](const uint64_t& i) -> const P_TYPE {return       -parametric[i].coeffs[2].dot(parametric[i].coeffs[3])  / polys[i].coeffs[0]; }
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

    void throw_arg_exception(const std::size_t ind, const ppl::vertex<P_TYPE>* const points) const {

        // out << "9-0-0870";
        std::string err = "control points are on top of each other!\nstarting from the "+std::to_string((ind * ppl::cubic)+1)+"'th control point:\n";
        for(std::size_t i{0}; i < ppl::cubic_points; ++i ){
            err+=std::to_string(i + (ind * ppl::cubic) +1);
            err += *(points+i);
        }

        std::cerr<< err << "\n";
        _invalid_argument("invalid argument");
    }


    void call_significant_fig_ascertain(void) {

        if(!(ppl::TOLERANCE>0) || ppl::TOLERANCE>=1)
            _invalid_argument("invalid tolerance value! the value of tolerance has to be in range  ] 1.0, 0.0 [ \n");

        int t_prec = -std::log10(ppl::TOLERANCE);
        if( ppl::prec_call(ppl::TOLERANCE) > std::numeric_limits< P_TYPE>::digits10 ){ 
            
            std::cerr.precision(t_prec+1); 
            std::cerr<< "the tolerance value of " << std::fixed << ppl::TOLERANCE << " is out of range for precision of type <data type " << ppl::__Tn<P_TYPE>() << "> \n";

            _logic_error("incompatible arguments");
        }
        
    }




    inline bool newton_mth(const ppl::sturmSeq<P_TYPE>& seq, P_TYPE val, 
                        const P_TYPE& a, const P_TYPE& b, 
                        ppl::real_roots<P_TYPE>& roots) const  
    {

        P_TYPE polyEvalu, derivEvalu;


        for (std::size_t i{0};;++i)
        {   
            polyEvalu = poly1d_solve_for(seq.poly[0], val); 

            if ( std::abs( polyEvalu ) <= TOLERZ){
                roots.push(val);
                return 1;
            }

            derivEvalu = poly1d_solve_for(seq.poly[1], val);

            if ( derivEvalu == 0.0 || i > NEWTON_THRES)  //       <<<<<<<<<<<<<<<<<<<    NEWTON'S METHOD FAILED!!
                return 0;  // >>>>  throw local maximum/minimum || interations overflow!
            

            val = val - ( polyEvalu / derivEvalu );


            if(val < a || val > b) //         <<<<<<<<<<<<<    NEWTON'S METHOD FAILED!!
                return 0;  //  >>>> throw wrong root
            
        }

    }


    inline uint8_t _signs_at
                (const ppl::sturmSeq<P_TYPE>& seq,
                const P_TYPE &val) const {

        uint16_t alters{0};
        bool priv, curr;
        priv = std::signbit(poly1d_solve_for(seq.poly[0],val));
    
        for (uint8_t i{1}; i < seq.len; ++i){
            curr =  std::signbit(poly1d_solve_for(seq.poly[i],val));
            alters+=(priv^curr);
            priv = curr;
        }

        return alters;
    }
    
    inline uint8_t roots_in_interval
                (const ppl::sturmSeq<P_TYPE>& seq,
                const P_TYPE& a, const P_TYPE& b) const {
        return std::abs(_signs_at(seq, a)-_signs_at(seq, b));
    }


    inline uint8_t __split(const ppl::sturmSeq<P_TYPE>& seq, 
                const P_TYPE& a, const P_TYPE& b, 
                const uint8_t& _rN, 
                ppl::real_roots<P_TYPE>& roots) const
    {
        if ((b - a) <= TOLERZ){
            roots.push( (a + b) / 2.0 ); 
            return _rN;
        }

        P_TYPE m_value = (a+b) / 2.0;

        if (_rN == 1)
        {

            P_TYPE l_evalu{poly1d_solve_for(seq.poly[0], a)}, r_evalu{poly1d_solve_for(seq.poly[0], b)};
            
            if ( l_evalu < 0.0 && r_evalu > 0.0)
            {

                /*
                    Newton's method is extremely fast to find a root, 
                    but if it FOR VERY RARE SITUATION failed to find a root in a certain number of iterations, 
                    then more likely it'll not find a root at all, or it could tend toward a wrong root!
                    For some situations such oscillating sequence it's a must to change 
                    the initial value by shrinking the interval using Bisection method.
                    
                */

               if(newton_mth(seq, m_value, a, b, roots))
                    return 1;

                else{
                    if (ppl::__sign(poly1d_solve_for(seq.poly[0], m_value)) == ppl::__sign(r_evalu))
                        return __split(seq, a, m_value, _rN, roots);
                    else
                        return __split(seq, m_value + TOLERZ, b, _rN, roots);
                }

            }else return 1;


        }else {

            uint8_t rootsN = roots_in_interval(seq, m_value + TOLERZ, b );

            if( rootsN != 0 ){
                rootsN = __split(seq, m_value + TOLERZ, b, rootsN, roots);
            }

            if(rootsN != _rN) 
                return rootsN + __split(seq, a, m_value, _rN -rootsN, roots);
            else return rootsN;
        }
    }


    inline bool poly_long_division
            (const ppl::poly1d<P_TYPE>& numer,
            const ppl::poly1d<P_TYPE>& denom, 
            ppl::poly1d<P_TYPE>& rem) const
    {

        memcpy(rem.coeffs, numer.coeffs, sizeof(P_TYPE) * ppl::quintic_Coeffs);
        rem.d = numer.d;
        P_TYPE q;

        uint8_t i, j;
        for ( ;denom.d<=rem.d; )
        {
            q = rem.coeffs[ppl::quintic-rem.d]/denom.coeffs[ppl::quintic-denom.d];
            rem.coeffs[ppl::quintic - rem.d ] = 0.0;
            i = ppl::quintic - rem.d + 1;
            j = i + denom.d;
            for (; i < j; ++i)
            {
                rem.coeffs[i] = rem.coeffs[i] - (q * denom.coeffs[rem.d - denom.d + i ] );
            }
            --rem.d;

            for (i = ppl::quintic- rem.d ; i < ppl::quintic_Coeffs; ++i){
                if (rem.coeffs[i] == 0.0){
                    if(rem.d == 0)
                        return 0;
                    else --rem.d;
                }
                    
                else
                    break;
            }
        }
        std::for_each(rem.coeffs+(ppl::quintic-rem.d),
                    rem.coeffs+ppl::quintic_Coeffs, 
                    [](P_TYPE& coeff){ coeff = -coeff;});

        return 1;
    }


    inline void construct_sturmPolys(ppl::sturmSeq<P_TYPE>& seq) const {
        
        std::size_t i;
        for (i = 1; i < ppl::quintic; ++i) // first derivative
            seq.poly[1].coeffs[i + 1] = (seq.poly[0].d - i) * seq.poly[0].coeffs[i];

        for (i = 1; seq.poly[i].d > 0; ++i){   
            if(poly_long_division(seq.poly[i-1], seq.poly[i], seq.poly[i+1]))
                seq.len++;
        }
    }

    inline P_TYPE poly1d_solve_for
            (const ppl::poly1d<P_TYPE>& poly,
            const P_TYPE& val) const
    {

        P_TYPE result = poly.coeffs[ppl::quintic-poly.d];

        std::for_each(poly.coeffs+(ppl::quintic-poly.d + 1), 
                        poly.coeffs+ppl::quintic_Coeffs,
                            [&](const P_TYPE& coeff){ result = result*val + coeff;}
                        );

        return result;
    }

    inline ppl::vertex<P_TYPE> poly3d_solve_for
            (const ppl::poly3d<P_TYPE>& poly, 
            const P_TYPE& val) const{

        return val*(val*(val*poly.coeffs[0]+poly.coeffs[1])+poly.coeffs[2]) + poly.coeffs[3];
    }

    inline void _call_projection(const ppl::vertex<P_TYPE> * const p)
    {
        
        _assert__(poly_num>0, "closest point was called on empty data! did you forget to load your data?\n");
        ppl::vertex<P_TYPE> rootsVertex;
        uint8_t _rN;

        curr_dist = (*p).sqr_dist( splines[0][0] );
        if (min_dist > curr_dist)
        {
            point_projection.closest = splines[0][0];
            point_projection.index = point_projection.parameter = 0;
            min_dist = curr_dist;
        }

        for (std::size_t i{0}; i < poly_num; ++i)
        {

            curr_dist = (*p).sqr_dist(splines[i][ppl::cubic]);
            if (min_dist > curr_dist)
            {
                point_projection.closest = splines[i][ppl::cubic];
                point_projection.index = i;
                point_projection.parameter = 1;
                min_dist = curr_dist;
            }

            ppl::sturmSeq<P_TYPE> seq;
            for(uint8_t j{1}; j<ppl::quintic_Coeffs; ++j)
                if(j<ppl::cubic)
                    seq.poly[0].coeffs[j] = polys[i].coeffs[j];
                else
                    seq.poly[0].coeffs[j] = polys[i].coeffs[j] 
                    + (deriv[i].coeffs[j-ppl::cubic].dot(*p) / polys[i].coeffs[0]);

            construct_sturmPolys(seq);


            _rN = roots_in_interval(seq, lowerB, upperB);

            if (_rN != 0){
                ppl::real_roots<P_TYPE> roots;
                __split(seq, lowerB, upperB, _rN, roots);

                for (std::size_t j{0}; j < roots.num; ++j){
                    rootsVertex = poly3d_solve_for(parametric[i], roots.zeros[j]);
                    curr_dist = (*p).sqr_dist(rootsVertex);

                    if (min_dist > curr_dist){
                        point_projection.index = i;
                        point_projection.parameter = roots.zeros[j];
                        min_dist = curr_dist;
                    }
                }
                // roots.clear();
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
    inline void __freem(T* &_alloc){
        if(_alloc != nullptr){
            delete[] _alloc;
            _alloc = nullptr;
        }
    }

public:    
    // cubic_path() = default;
    cubic_path(): 
            splines{nullptr}, control_points{nullptr},   
            polys{nullptr},parametric{nullptr}, 
            deriv{nullptr},
            poly_num{0}, points_num{0},
            min_dist{std::numeric_limits<P_TYPE>::max()},
            curr_dist{min_dist} {
                _assert__(std::is_floating_point<P_TYPE>::value, "instantiation of ppl::cubic_path can only be with floating-point types!\n");}

    cubic_path(const ppl::vertex<P_TYPE>* const points, const uint64_t& _size){
        _assert__(std::is_floating_point<P_TYPE>::value, "instantiation of ppl::cubic_path can only be with floating-point types!\n");
        routing(points, _size);
    }
       
    virtual ~cubic_path() { cleanUp(); }


#ifdef __CONCURRENCY__

    void closest_point(const ppl::vertex<P_TYPE> * const p,  ppl::projection<P_TYPE> * const projection_ptr)
    {
        _call_projection(p);
        if(point_projection.parameter != 1 && point_projection.parameter != 0){
            projection_ptr->closest = poly3d_solve_for(parametric[point_projection.index], point_projection.parameter);
        }else{
            projection_ptr->closest = point_projection.closest;
        }
        
        projection_ptr->index = point_projection.index;
        projection_ptr->parameter = point_projection.parameter;
        projection_ptr->dist = projection_ptr->closest.dist(*p);

        min_dist = std::numeric_limits<P_TYPE>::max();
    }

#endif

    void routing(const ppl::vertex<P_TYPE>* const points, const uint64_t& _size)
    {      
        _assert__( (_size -1)%ppl::cubic == 0 && _size > ppl::cubic , "incompatible number of control points!");
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
            memcpy(splines[j], points+i, sizeof(ppl::vertex<P_TYPE>)*ppl::cubic_points);
        }

        polys = new ppl::poly1d<P_TYPE>[poly_num];

        parametric = new ppl::poly3d<P_TYPE>[poly_num];
        deriv = new ppl::deriv3d<P_TYPE>[poly_num];

        for(i=0, j=0; i < term; ++j, i+=ppl::cubic){
            extract_poly(j, points+i);
        }
        
        min_dist = std::numeric_limits<P_TYPE>::max();
    }

    ppl::projection<P_TYPE> closest_point(const ppl::vertex<P_TYPE> * const p)

    {
        _call_projection(p);
        min_dist = std::numeric_limits<P_TYPE>::max();
        if(point_projection.parameter != 1 && point_projection.parameter != 0){
            point_projection.closest = poly3d_solve_for(parametric[point_projection.index], point_projection.parameter);
        }

        point_projection.dist = point_projection.closest.dist(*p);
        return point_projection;
    
    }

};
#ifdef __CONCURRENCY__

template<typename  P_TYPE>  struct thrStr{

    static void* _task(void* argv){
        
        ppl::cubic_path<P_TYPE>* path = (ppl::cubic_path<P_TYPE>*) ((void**)argv)[0];

        path->closest_point(  (ppl::vertex<P_TYPE> *) ((void**)argv)[1], (ppl::projection<P_TYPE> *) ((void**)argv)[2] );

        pthread_exit(nullptr);
    }
};

#endif
template< typename P_TYPE> class point_projection

{   
    static_assert(std::is_floating_point<P_TYPE>::value, "instantiation of ppl::point_projection can only be with floating-point types!\n");

#ifdef __EXTERNAL_TRACK_LOADING__

    ppl::vertex<P_TYPE>* points{nullptr};
    uint64_t _size{0};

    void* addr{nullptr};
    std::istringstream s_stream;
    signed f_descriptor;
    struct stat buf;


#endif

#ifdef __CONCURRENCY__

    std::size_t _thrN{0};
    std::vector<std::tuple<uint64_t, uint64_t>> jobs_intervals;
    std::vector<ppl::projection<P_TYPE>> thr_verts;

    ppl::cubic_path<P_TYPE>* _track_strips{nullptr};

    uint64_t splinesN;
    _channel** thrd{nullptr};
    
#else

    ppl::cubic_path<P_TYPE> _track;
    
#endif

#ifdef __CONCURRENCY__
    void build_intervals(const unsigned& tks_per_thr, const unsigned& __size)
    {
        jobs_intervals.clear();
        std::tuple<std::size_t, std::size_t> tup;

        std::get<0>(tup) = 0;
        std::get<1>(tup) = (tks_per_thr * ppl::cubic) + 1;
        jobs_intervals.push_back(tup);

        for (unsigned j{1}; j < _thrN; ++j){

            std::get<0>(tup) = std::get<0>(jobs_intervals[j-1]) + std::get<1>(jobs_intervals[j-1]) - 1;
            std::get<1>(tup) =  (tks_per_thr * ppl::cubic) + 1;
            jobs_intervals.push_back(tup);
        }

        std::get<1>(jobs_intervals[jobs_intervals.size() - 1]) = __size - std::get<0>(jobs_intervals[jobs_intervals.size() - 1]);
    }
#endif

    template< typename T > inline void __freem(T* &_alloc){
        if(_alloc != nullptr){
            delete[] _alloc;
            _alloc = nullptr;
        }
    }

    void cleanUp(void){

#ifdef __CONCURRENCY__

        __freem(_track_strips);
        __freem(thrd);
#endif

#ifdef __EXTERNAL_TRACK_LOADING__
        
        __freem(points);
#endif

    }

public:
    point_projection() = default;
    virtual ~point_projection(){
        
        cleanUp();

    }


#ifdef __EXTERNAL_TRACK_LOADING__

    uint64_t num_of_p(char *addr, const uint64_t& len){

        std::stringstream s_stream;
        s_stream.rdbuf()->pubsetbuf(addr, len);
        return std::distance(std::istream_iterator<std::string>(s_stream), std::istream_iterator<std::string>());
    }


    inline static bool does_exist(const fs::path& _dir, fs::file_status _status = fs::file_status{})
    {
	    return (fs::status_known(_status) ? fs::exists(_status) : fs::exists(_dir));

    }

    void routing(const std::string& _dir){
        cleanUp();
#ifdef __linux__

        if(!does_exist(_dir))
            throw std::invalid_argument("cannot find file <"+_dir+">\n" );

        f_descriptor = open(_dir.c_str(), O_RDONLY, S_IRUSR|S_IWUSR);
        
        if( f_descriptor < 0)
            throw std::runtime_error("could not open file <"+_dir+">\n");

        if(fstat(f_descriptor,&buf) < 0){
            close(f_descriptor);
            throw std::runtime_error("not able to get file size <"+_dir+">\n");
        }

        addr = mmap(NULL, buf.st_size, PROT_READ, MAP_PRIVATE, f_descriptor, 0);

        if(addr == MAP_FAILED){
            close(f_descriptor);
            throw std::runtime_error("an exception occurred while mapping the file <"+_dir+">\n");
        }


        s_stream.rdbuf()->pubsetbuf ( reinterpret_cast<char*>(addr), buf.st_size );
        _size =  num_of_p((char*)addr, buf.st_size);

        if( _size%ppl::cubic != 0){
            munmap(addr, buf.st_size);
            close(f_descriptor);
            throw std::logic_error("incompatible size of read data in the file <"+_dir+">, "+ std::to_string(_size)+" were read!");
        }
        
        _size /= ppl::cubic;
        points = new ppl::vertex<P_TYPE>[_size];

        P_TYPE x{0}, y{0}, z{0};
        std::size_t i{0};   
        while(s_stream >> x >> y >> z)
            points[i++] = ppl::vertex<P_TYPE>{ x, y, z};
        s_stream.clear();
        
        munmap(addr, buf.st_size);
        close(f_descriptor);

        _assert__((_size -1)%ppl::cubic == 0 && _size > ppl::cubic, "incompatible number of control points!");


#if defined __CONCURRENCY__   

        drafting_concur_attrib(this->points, this->_size);

#else
        _track.routing(this->points, this->_size);
         
#endif

#else
        
#error "loading control points from external file is not supported on this system!"

#endif

    }
    
    point_projection(const std::string& _dir){
        routing(_dir);
    }

#endif

    point_projection(const ppl::vertex<P_TYPE>* const _points,const uint64_t& __size) {
        
        routing(_points, __size);
    }

    void routing(const ppl::vertex<P_TYPE>* const _points, const uint64_t& __size){
        _assert__( (__size -1)%ppl::cubic == 0 && __size > ppl::cubic, "incompatible number of control points!" );

        cleanUp();

#if defined __CONCURRENCY__ 

        drafting_concur_attrib(_points, __size);

#else
        _track.routing(_points, __size);
         
#endif
        
    }


#if defined __CONCURRENCY__ 
    
    void drafting_concur_attrib(const ppl::vertex<P_TYPE>* const _points, const uint64_t& __size){

        splinesN = (__size-1)/ppl::cubic;
      

#if defined(_WIN32) || defined(WIN32) 

        _SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        _thrN = sysinfo.dwNumberOfProcessors;

#elif defined __linux__

        _thrN = sysconf(_SC_NPROCESSORS_CONF); 

#else
        
#error "multithreading support is not supported on this system!"

#endif
        _track_strips = new ppl::cubic_path<P_TYPE>[_thrN];

        if(_thrN>splinesN) 
            _thrN = splinesN;
        thrd = new _channel*[_thrN];

        thr_verts.resize(_thrN);

        build_intervals(splinesN/_thrN, __size);
        
        for(std::size_t i{0}; i < jobs_intervals.size(); ++i)
        {
            _track_strips[i].routing(_points+std::get<0>(jobs_intervals[i]), std::get<1>(jobs_intervals[i]));
        }
    }

#endif


    ppl::projection<P_TYPE> localize(const ppl::vertex<P_TYPE> * const p)
    
    {

#if defined __CONCURRENCY__

        int32_t i, min_ind, _stride{0};

        for (i = _thrN-1; i >=0 ; --i)
            thrd[i] = new _channel( ppl::thrStr<P_TYPE>::_task, 3, &_track_strips[i], p, &thr_verts[i]);

        for (i = _thrN-1; i >=0 ; --i)
            thrd[i]->join(nullptr);
        for (i = _thrN-1; i >=0 ; --i)
            delete thrd[i];

        min_ind = std::min_element(thr_verts.begin(), thr_verts.end(),
                        [](const ppl::projection<P_TYPE> v1, const ppl::projection<P_TYPE> v2)
                        ->bool{return v1.dist<v2.dist;} ) - thr_verts.begin(); 


        if(min_ind == 0)
            return thr_verts[0];
            
        for( i = 0; i < min_ind; ++i)
            _stride+= ( std::get<1>(jobs_intervals[i]) -1)/ppl::cubic;

        thr_verts[min_ind].index += _stride;

        return  thr_verts[min_ind];
#else 

        return _track.closest_point(p);
        
#endif
    }
};



} // namespace ppl





#endif   //  QUINTIC_POLY_HPP

