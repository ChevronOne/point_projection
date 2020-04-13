
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



template<typename P_TYPE> PPL_FUNC_DECL ppl::LERPer::__state<P_TYPE> 
state_call(const std::vector<ppl::vertex<P_TYPE>>& data, 
        const ppl::vertex<P_TYPE>* const controlP,
        const ppl::LERPer::_sector& _sec)
{
    ppl::cubic_path<P_TYPE> path(controlP, ppl::cubic_points);
    ppl::LERPer::__state<P_TYPE> STATE;
    for(std::size_t i{_sec.f}; i<=_sec.l; ++i){
        P_TYPE dist = path.closest_point(&data[i]).dist;
        if (dist > STATE.ERR){
            STATE.ERR = dist;
            STATE.IND = i;
        }
    }
    return STATE;
}


template<typename P_TYPE> ppl::LERPer::__state<P_TYPE>
extractB_path(const std::vector<ppl::vertex<P_TYPE>>& data, 
        std::vector<ppl::vertex<P_TYPE>>& controlP,
        const ppl::LD& TOLEZ)
{
    if(data.size()<ppl::cubic_points)
        ppl_invalid_argument("no enough data to fit!");
    ppl::LERPer::call_significant_fig_ascertain<P_TYPE>(TOLEZ);
    std::vector<uint64_t> strips;
    ppl::LERPer::_winnow(data, strips);

    uint64_t actually_size = controlP.size();
    controlP.reserve(controlP.size()+strips.size());
    controlP.emplace_back(data.front());

    ppl::LERPer::__state<P_TYPE> outcome, _clarific;
    boost::progress_display prog_bar(data.size(), std::cout, "adjusting curves: \n");


    uint64_t it{strips.size() <= ppl::cubic_Coeffs? strips.size()-1 : ppl::cubic};
    ppl::LERPer::_sector _sec{strips[0], strips[it], strips[it]+1};
    bool STRIP_SUCC{0};
    for(;_sec.f!=_sec.l;){ 

        STRIP_SUCC = ppl::LERPer::sector_tuneUP(data, controlP, outcome, _sec, TOLEZ);

        if(STRIP_SUCC){

            for(;strips[it-1]>_sec.l; --it){}
                
            it = it+ppl::cubic<strips.size()? it+ppl::cubic:strips.size()-1;
            prog_bar = _sec.l;

            _sec.f = _sec.l;
            _sec.l = strips[it];

            if(outcome.ERR >_clarific.ERR)
                _clarific = outcome;
        }
        else _sec.l = outcome.IND; 
        _sec.len = _sec.l-_sec.f+1; 
    }

    prog_bar = data.size();
    std::cout << "\n data size: " << data.size() 
              << "\n number of fitted curves: "
              << (controlP.size()-actually_size-1)/ppl::cubic
              << "\n max error: " << _clarific.ERR
              << "\n"; 

    return _clarific;
}        



template<typename P_TYPE> ppl::LERPer::__state<P_TYPE>
fitSingle(const std::vector<ppl::vertex<P_TYPE>>& data, 
        std::vector<ppl::vertex<P_TYPE>>& controlP,
        const P_TYPE& TOLEZ) 
{
    if(data.size()<ppl::cubic_points)
        ppl_invalid_argument("no enough data to fit!");
    ppl::LERPer::call_significant_fig_ascertain<P_TYPE>(TOLEZ);
    ppl::LERPer::__state<P_TYPE> _clarific = ppl::LERPer::_attempt_to_fit(data, controlP, 
                    {0, data.size()-1, data.size()}, 
                    TOLEZ);

    std::cout << "\n data size: " << data.size() 
              << "\n max error: " << _clarific.ERR
              << "\n";
    return _clarific;
}



template<typename P_TYPE> PPL_FUNC_DECL ppl::LERPer::__state<P_TYPE> 
_attempt_to_fit(const std::vector<ppl::vertex<P_TYPE>>& data, 
        std::vector<ppl::vertex<P_TYPE>>& controlP,
        const ppl::LERPer::_sector& _sec,
        const ppl::LD& TOLEZ)
{
    if(_sec.len<ppl::cubic)
        ppl_invalid_argument("no enough data to fit!");

    controlP.resize(ppl::cubic_points);
    ppl::vertex<P_TYPE> fitted_CP[ppl::cubic_points];
    fitted_CP[0]= controlP[0] = data[_sec.f];
    fitted_CP[3]= controlP[3] = data[_sec.l]; 


    ppl::vertex<P_TYPE> V1{left_tan(data, _sec.f)}, V2{right_tan(data, _sec.l)};
    std::vector<P_TYPE> para(_sec.len);
    ppl::LERPer::_centripetal_par(data, para, _sec);
    ppl::cubic_path<P_TYPE> path;
    ppl::LERPer::__state<P_TYPE> STATE, OUTCOME;
    double CONV{std::numeric_limits<double>::max()};
    P_TYPE P_DIST{std::numeric_limits<P_TYPE>::max()};


    uint64_t k{0};
    for(k = 0; CONV >= VERGENCE && P_DIST >= TOLEZ && k<ppl::IT_OVR_FLW; ++k)
    {

        P_TYPE C11{0}, C12{0}, C22{0}, X1{0}, X2{0};

        uint64_t i{0},j{_sec.f};
        for (i = 0; i < _sec.len; ++j, ++i)
        {
            ppl::vertex<P_TYPE> A1 {basisF<P_TYPE>[1](para[i]) * V1},
                                A2 {basisF<P_TYPE>[2](para[i]) * V2};

            C11 += A1.dot(A1);
            C12 += A1.dot(A2);
            C22 += A2.dot(A2);

            ppl::vertex<P_TYPE> V0 = data[j]
                     - ((basisF<P_TYPE>[0](para[i]) + basisF<P_TYPE>[1](para[i])) * fitted_CP[0]
                     +  (basisF<P_TYPE>[2](para[i]) + basisF<P_TYPE>[3](para[i])) * fitted_CP[3]);

            X1 += A1.dot(V0);
            X2 += A2.dot(V0);
        }

        P_TYPE det{(C11 * C22 - C12 * C12)}; 
        P_TYPE  alph1{ static_cast<P_TYPE>((det==0.0)? 0.0 : (C22 * X1 - C12 * X2) / det) }, 
                alph2{ static_cast<P_TYPE>((det==0.0)? 0.0 : (C11 * X2 - C12 * X1) / det) }; 

        if(alph1<=0.0||alph2<=0.0){
            alph1 = fitted_CP[0].dist(fitted_CP[3])/2.0;
            fitted_CP[1] = alph1 * V1 + fitted_CP[0];
            fitted_CP[2] = alph1 * V2 + fitted_CP[3];

        }else{
            fitted_CP[1] = alph1 * V1 + fitted_CP[0];
            fitted_CP[2] = alph2 * V2 + fitted_CP[3];
        }

        path.routing(fitted_CP, ppl::cubic_points);
        STATE.ERR = 0.0;
        for (j = _sec.f , i = 0; j <= _sec.l; ++j, ++i){
            ppl::projection<P_TYPE> C_STATE = path.closest_point(&data[j]);
            para[i] = C_STATE.parameter;
            if (C_STATE.dist > STATE.ERR){
                STATE.ERR = C_STATE.dist;
                STATE.IND = j;
            }
        }

        CONV = P_DIST - STATE.ERR;
        if(STATE.ERR<P_DIST){
            controlP[1] = fitted_CP[1];
            controlP[2] = fitted_CP[2];
            OUTCOME = STATE;
        }
        P_DIST = STATE.ERR;
    }

    return OUTCOME;
}



template<typename P_TYPE> PPL_FUNC_DECL void
_centripetal_par(const std::vector<ppl::vertex<P_TYPE>>& data, 
        std::vector<P_TYPE>& para,
        const ppl::LERPer::_sector& _sec) 
{
    if(data.size()<2)
        ppl_invalid_argument("no data to generate parameter values!");

    para[0]=0;
    for(uint64_t i{_sec.f+1}, j=1; i<=_sec.l; ++i, ++j)
        para[j] = para[j-1] + data[i].dist(data[i-1]);

    for(P_TYPE& p:para)
        p = p/para.back();
}



template<typename P_TYPE> PPL_FUNC_DECL ppl::vertex<P_TYPE> 
right_tan(const std::vector<ppl::vertex<P_TYPE>>& data, const uint64_t& ind)
{
    ppl_assert__(ind>0, "right tangent called on the most left data point!");
    ppl_assert__(ind<data.size(), "out of range index of right tangent!");

    if(ind == (data.size()-1))
        return ppl::normalize(*std::next(data.rbegin()) - *data.rbegin());
    
    else
        return ppl::normalize(data[ind-1] - data[ind+1]);
    
}



template<typename P_TYPE> PPL_FUNC_DECL bool 
sector_tuneUP(const std::vector<ppl::vertex<P_TYPE>>& data, 
        std::vector<ppl::vertex<P_TYPE>>& controlP,
        ppl::LERPer::__state<P_TYPE>& outcome,
        const ppl::LERPer::_sector& _sec,
        const ppl::LD& TOLEZ)
{
    ppl_assert__(_sec.len>1, "no data to fit!\n");

    if(_sec.len == 2){

        P_TYPE DIST{data[_sec.f].dist(data[_sec.l])/ static_cast<P_TYPE>(3.0)};

        controlP.emplace_back(left_tan(data, _sec.f)*DIST);
        controlP.emplace_back(right_tan(data, _sec.l)*DIST);
        controlP.emplace_back(data[_sec.l]);

        outcome = {_sec.f, 0.0}; 

        return 1;

    }else{

        std::vector<ppl::vertex<P_TYPE>> p_container;
        outcome = ppl::LERPer::_attempt_to_fit(data, p_container, 
                        _sec, TOLEZ);

        if(outcome.ERR<=TOLEZ){

            controlP.insert(controlP.end(), std::next(p_container.begin()), p_container.end());
            return 1;

        }else return 0;

    }

}


template<typename P_TYPE> PPL_FUNC_DECL ppl::vertex<P_TYPE> 
left_tan(const std::vector<ppl::vertex<P_TYPE>>& data, const uint64_t& ind)
{
    ppl_assert__(ind<(data.size()-1), "left tangent called on the most right data point!");
    ppl_assert__(ind>=0, "out of range index of left tangent!");

    if(ind == 0)
        return ppl::normalize(*std::next(data.begin()) - *data.begin());

    else
        return ppl::normalize(data[ind+1] - data[ind-1]);

}



template<typename P_TYPE> PPL_FUNC_DECL
void _winnow(const std::vector<ppl::vertex<P_TYPE>>& data,
             std::vector<uint64_t>& strips)
{
    if(data.size()<=ppl::cubic_points)
        ppl_invalid_argument("no enough data to split!");
    strips.clear();
    strips.push_back(0);

    bool dir_ch{0};
    for(uint64_t i{1}; i < data.size()-1; ++i){

        if( (data[i-1].x > data[i].x)^(data[i].x > data[i+1].x) ||

            (data[i-1].y > data[i].y)^(data[i].y > data[i+1].y) ||

            (data[i-1].z > data[i].z)^(data[i].z > data[i+1].z) )
            {dir_ch^=1;}

        if(dir_ch){
            strips.push_back(i);
            dir_ch^=1;
        } 
    }

    if(strips.back() != data.size()-1)
        strips.push_back(data.size()-1); 
}



template<typename P_TYPE> PPL_FUNC_DECL void 
call_significant_fig_ascertain(const ppl::LD& TOLEZ){
        if(TOLEZ<=0)
                ppl_invalid_argument("invalid tolerance value! the value of tolerance has to be greater than 0.0 \n");

		std::streamsize t_prec = static_cast<std::streamsize>( -std::log10(TOLEZ) );
        if(ppl::prec_call(TOLEZ) > std::numeric_limits< P_TYPE>::digits10 ){

                std::streamsize def_prec = std::cerr.precision();
                std::cerr.precision(t_prec +1);
                std::cerr<< "the tolerance value of " << std::fixed << TOLEZ << " is out of range for precision of type <data type " <<  ppl::__Tn<P_TYPE>()  << "> \n";
                std::cerr.precision(def_prec);
                ppl_logic_error("incompatible arguments");
        }
}







