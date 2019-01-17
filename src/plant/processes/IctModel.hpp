/**
 * @file ecomeristem/plant/thermal-time/Model.hpp
 * @author The Ecomeristem Development Team
 * See the AUTHORS or Authors.txt file
 */

/*
 * Copyright (C) 2005-2017 Cirad http://www.cirad.fr
 * Copyright (C) 2012-2017 ULCO http://www.univ-littoral.fr
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef ICT_MODEL_HPP
#define ICT_MODE_HPP

#include <defines.hpp>

namespace model {

class IctModel : public AtomicModel < IctModel >
{
public:    
    enum internals { IS_COMPUTED, SURVIVED, IC_1, IC, TOTAL, N };
    enum externals { BOOL_CROSSED_PHYLLO, IC_plant };


    IctModel() {
        //    computed variables
        Internal(IS_COMPUTED, &IctModel::_is_computed);
        Internal(SURVIVED, &IctModel::_survived);
        Internal(IC_1, &IctModel::_ic_1);
        Internal(IC, &IctModel::_ic);
        Internal(TOTAL, &IctModel::total);
        Internal(N, &IctModel::n);

        //    external variables
        External(BOOL_CROSSED_PHYLLO, &IctModel::_bool_crossed_phyllo);
        External(IC_plant, &IctModel::_ic_plant);
    }

    virtual ~IctModel()
    {}

    void compute(double t, bool /* update */) {
        if (t != _parameters.beginDate) {
            if(_is_computed) {
                _ic_ = _ic_;
                _ic = _ic;
                _survived = _survived;
                _is_computed = _is_computed;
                return;
            } else {
                if(_bool_crossed_phyllo >= 0) {
                    //test si on va survivre
                    if(_ic > Ict) {
                        _survived = true;
                    } else {
                        _survived = false;
                    }
                    _is_computed = true;
                } else {
                    //declaration variables
                    total = 0;
                    n = 0;
                    double mean = 0;
                    _ic_1 = _ic;

                    //ajout des valeurs du jour
                    _ic_.push_back(_ic_plant);

                    //calcul ic
                    for(int i = 0; i < _ic_.size(); i++) {
                        total += _ic_[i];
                        ++n;
                    }

                    if(n != 0) {
                        mean = total/n;
                    } else {
                        mean = _ic_1;
                    }

                    double tmp = std::min(5., mean);
                    if(tmp == 0) {
                        _ic = 0.0001;
                    } else {
                        _ic = tmp;
                    }
                }
            }
        }
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        //parameters
        _parameters = parameters;
        Ict = _parameters.get("Ict");

        //internals
        _is_computed = false;
        _survived = true;
        total = 0;
        n = 1;
        _ic_1 = 0.0001;
        _ic = 0.0001;
    }

private:
    ecomeristem::ModelParameters _parameters;
    //parameters
    double Ict;

    //internals
    bool _is_computed;
    bool _survived;

    double _ic_1;
    double _ic;
    double total;
    double n;

    std::vector<double> _ic_;

    //externals
    double _bool_crossed_phyllo;
    double _ic_plant;


};

} // namespace model
#endif
