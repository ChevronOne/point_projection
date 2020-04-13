

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


#ifndef PPL_QUINTIC_POLY_HPP
#define PPL_QUINTIC_POLY_HPP


#if defined(__CONCURRENCY__) &&  !defined(_WIN32) && !defined(WIN32) && !defined(__linux__)
#error "multithreading support is not supported on this system!"
#elif defined(__EXTERNAL_TRACK_LOADING__) && !defined(__linux__)

#error "loading control points from external file is not supported on this system!"
#endif





namespace ppl{

#ifdef __CONCURRENCY__

template<typename  P_TYPE>  struct thrStr{

    static void* _task(void* argv){
        
        ppl::cubic_path<P_TYPE>* path = (ppl::cubic_path<P_TYPE>*) ((void**)argv)[0];

        path->closest_point(  (ppl::vertex<P_TYPE> *) ((void**)argv)[1],
                     (ppl::projection<P_TYPE> *) ((void**)argv)[2] );

        pthread_exit(nullptr);
        return nullptr;
    }
};

#endif

template< typename P_TYPE> class point_projection

{   
	static_assert(std::numeric_limits<P_TYPE>::is_iec559,
		"instantiation of ppl::point_projection can only be with floating-point types!\n");

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
    void build_intervals(const unsigned& tks_per_thr, 
        const unsigned& __size)
    {
        jobs_intervals.clear();
        std::tuple<std::size_t, std::size_t> tup;

        std::get<0>(tup) = 0;
        std::get<1>(tup) = (tks_per_thr * ppl::cubic) + 1;
        jobs_intervals.push_back(tup);

        for (unsigned j{1}; j < _thrN; ++j){

            std::get<0>(tup) = std::get<0>(jobs_intervals[j-1])
                    + std::get<1>(jobs_intervals[j-1]) - 1;
            std::get<1>(tup) =  (tks_per_thr * ppl::cubic) + 1;
            jobs_intervals.push_back(tup);
        }

        std::get<1>(jobs_intervals[jobs_intervals.size() - 1]) = __size 
                - std::get<0>(jobs_intervals[jobs_intervals.size() - 1]);
    }
#endif

    template< typename T > PPL_FUNC_DECL void __freem(T* &_alloc){
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
        return std::distance(std::istream_iterator<std::string>(s_stream), 
                std::istream_iterator<std::string>());
    }


    PPL_FUNC_DECL static bool 
    does_exist(const fs::path& _dir, 
            fs::file_status _status = fs::file_status{})
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

        addr = mmap(NULL, buf.st_size, 
                PROT_READ, MAP_PRIVATE, f_descriptor, 0);

        if(addr == MAP_FAILED){
            close(f_descriptor);
            throw std::runtime_error("an exception occurred while mapping the file <"+_dir+">\n");
        }

        s_stream.rdbuf()->pubsetbuf ( reinterpret_cast<char*>(addr), 
                                        buf.st_size );
        _size =  num_of_p((char*)addr, buf.st_size);

        if( _size%ppl::cubic != 0){
            munmap(addr, buf.st_size);
            close(f_descriptor);
            throw std::logic_error("incompatible size of read data in the file <"+_dir+">, "
                    + std::to_string(_size)+" were read!");
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

        ppl_assert__((_size -1)%ppl::cubic == 0 && _size > ppl::cubic, 
            "incompatible number of control points!");


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

    point_projection(ppl::vertex<P_TYPE> const * const _points,
        const uint64_t& __size) {
        
        routing(_points, __size);
    }

    void routing(ppl::vertex<P_TYPE> const * const _points, 
                const uint64_t& __size){
        ppl_assert__( (__size -1)%ppl::cubic == 0 && __size > ppl::cubic, 
                "incompatible number of control points!" );

        cleanUp();

#if defined __CONCURRENCY__ 

        drafting_concur_attrib(_points, __size);

#else
        _track.routing(_points, __size);
         
#endif
        
    }


#if defined __CONCURRENCY__ 
    
    void drafting_concur_attrib(ppl::vertex<P_TYPE> const * const _points, 
        const uint64_t& __size){

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
            _track_strips[i].routing(_points+std::get<0>(jobs_intervals[i]), 
                        std::get<1>(jobs_intervals[i]));
        }
    }

#endif


    ppl::projection<P_TYPE> localize(ppl::vertex<P_TYPE> const * const p)
    
    {

#if defined __CONCURRENCY__

        int32_t i, min_ind, _stride{0};

        for (i = _thrN-1; i >=0 ; --i)
            thrd[i] = new _channel( ppl::thrStr<P_TYPE>::_task, 3, 
                        &_track_strips[i], p, &thr_verts[i]);

        for (i = _thrN-1; i >=0 ; --i)
            thrd[i]->join(nullptr);

        for (i = _thrN-1; i >=0 ; --i)
            delete thrd[i];

        min_ind = std::min_element(thr_verts.begin(), 
                        thr_verts.end(),
                        [](const ppl::projection<P_TYPE> v1, 
                                const ppl::projection<P_TYPE> v2)
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


#endif   //  PPL_QUINTIC_POLY_HPP



