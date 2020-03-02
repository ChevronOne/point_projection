


#ifndef PPL_PTHRS_LIB_HPP
#define PPL_PTHRS_LIB_HPP


#include <cstddef>
#include <stdarg.h>
#include <pthread.h>
#include <errno.h>
#include <memory>





namespace ppl
{

class _channel
{
    pthread_t id;
    pthread_attr_t attrib;

    std::unique_ptr<void*[]> argv;


    va_list ptrVar_list;

    bool joinable{1};

public:
    
    _channel() = default;
    _channel(void *(*_task)(void *), std::size_t argc, ...) 
    {

        argv = std::make_unique<void*[]> (argc);

        va_start(ptrVar_list, argc);

        for(std::size_t i = 0; i < argc; ++i)
            argv[i] = va_arg( ptrVar_list, void * );


        va_end(ptrVar_list);

        pthread_attr_init(&attrib);
        pthread_create(&id, &attrib, _task, argv.get());
        
    }
    int join(void **_returns)
    {
        if(joinable){
            joinable ^= 1;
            return pthread_join(id, _returns);
        }

        return EINVAL;
    }


    int detach()
    {
        if(joinable){
            joinable ^= 1;
            return pthread_detach(id);
        }

        return EINVAL;
    }

    bool isJoinable(){return joinable;}
    ~_channel(){}

};

} // namespace ppl







#endif // PPL_PTHRS_LIB_HPP


