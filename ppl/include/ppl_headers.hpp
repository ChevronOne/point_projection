

#ifndef PPL_HEADERS_HPP
#define PPL_HEADERS_HPP



   #include <cstddef>
    #include <iostream>
    #include <iomanip>
    #include <string.h>
    #include <limits>
    #include <cmath>
    #include <string>
    #include <typeinfo>
    #include <type_traits>
    #include <algorithm>
    #include <functional>
    
    #include <stdio.h>
    #include <errno.h>

namespace ppl
{

    #define CONST const
    
    typedef std::size_t UNS;
    typedef long double LD;


    CONST UNS quintic{5};
    CONST UNS quartic{4};
    CONST UNS quintic_Coeffs{6};
    CONST UNS cubic_points{4};
    CONST UNS cubic{3};
    CONST UNS cubic_Coeffs{4};
    CONST UNS quadratic{2};
    CONST UNS quadratic_Coeffs{3};
    CONST UNS IT_OVR_FLW{3'000};
    

    CONST LD TOLERANCE(0.00001L);   //  user defined tolerance as an allowed margin of error for closest point 
 

    #define VERGENCE 0.001




    #define NUM_OF(arr)  sizeof(arr) / sizeof(*arr)

} // namespace ppl


    #ifdef __EXTERNAL_TRACK_LOADING__
        // #include <utility>
        #include <fstream>

#if __has_include(<filesystem>)
        #include <filesystem>
        namespace fs = std::filesystem;

#elif __has_include(<experimental/filesystem>)
        #include <experimental/filesystem>
        namespace fs = std::experimental::filesystem::v1;

#else
#error "the under use compiler dose not support C++ filesystem!"

#endif
        // #include <experimental/string_view>
        #include <sstream>


        #include <sys/mman.h>
        #include <unistd.h>
        #include <sys/stat.h>
        #include <stdlib.h>
        #include <fcntl.h>
        #include <iterator>





    #endif

    #ifdef __CONCURRENCY__

#if __has_include(<unistd.h>)
        #include <unistd.h>
#endif
        // #include <numeric>
        #include <vector>
        #include "ppl_pthlib.hpp"
        #include <cstdlib>

        #include <tuple>

#if __has_include(<minwindef.h>)
        #include <minwindef.h>
#endif

        #if defined(_WIN32) || defined(WIN32)
            #include<windows.h>
            #include<Windows.h>
        #endif

    #endif


#endif // PPL_HEADERS_HPP

