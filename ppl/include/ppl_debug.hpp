

#ifndef PPL_DEBUG_HPP
#define PPL_DEBUG_HPP

#include "ppl_headers.hpp"

namespace ppl
{


#if defined(_MSC_VER)

#define BUF_SIZE_ppl 0x01F4

#define error_msg []() -> const char* { char _BUF[BUF_SIZE_ppl]; strerror_s(_BUF, BUF_SIZE_ppl, errno); return _BUF; }

#define reset_errNo() (errno == 0 ? "None" : error_msg())

#else

#define reset_errNo() (errno == 0 ? "None" : strerror(errno))

#endif

#define log_err(M, ...) fprintf(stderr,\
        "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__,\
        reset_errNo(), ##__VA_ARGS__)

#define _assert__(A, M, ...) if(!(A)) {\
    log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _out_of_range(M, ...) {\
    log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _logic_error(M, ...) {\
    log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _invalid_argument(M, ...) {\
    log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _runtime_error(M, ...) {\
    log_err(M, ##__VA_ARGS__); errno=0; exit(EXIT_FAILURE); }

#define _incompatible_platform(M, ...) {\
    std::cerr << "[ERROR] " << "errno:"<< reset_errNo() << ",  "<< M << "\n"; errno=0; exit(EXIT_FAILURE); }


} // namespace ppl

#endif // PPL_DEBUG_HPP




