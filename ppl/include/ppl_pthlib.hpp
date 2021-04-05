

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



#ifndef PPL_PTHRS_LIB_HPP
#define PPL_PTHRS_LIB_HPP




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

    bool isJoinable(){return joinable;};
    ~_channel(){}

};

} // namespace ppl


#endif // PPL_PTHRS_LIB_HPP





