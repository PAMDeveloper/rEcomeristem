/**
 * @file ecomeristem/root/Model.hpp
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

#ifndef ROOT_MODEL_HPP
#define ROOT_MODEL_HPP

#include "defines.hpp"

namespace model {

class RootModel : public AtomicModel < RootModel >
{

public:
    enum internals { ROOT_DEMAND_COEF, ROOT_DEMAND, SURPLUS, LAST_ROOT_DEMAND, ROOT_BIOMASS, LAST_VALUE };

    enum externals { LEAF_DEMAND_SUM, LEAF_LAST_DEMAND_SUM, INTERNODE_DEMAND_SUM,
                     INTERNODE_LAST_DEMAND_SUM, PLANT_PHASE, PLANT_STATE,
                     CULM_SURPLUS_SUM };


    RootModel() {
        //    computed variables
        Internal(ROOT_DEMAND_COEF, &RootModel::_root_demand_coef);
        Internal(ROOT_DEMAND, &RootModel::_root_demand);
        Internal(SURPLUS, &RootModel::_surplus);
        Internal(LAST_ROOT_DEMAND, &RootModel::_last_root_demand);
        Internal(ROOT_BIOMASS, &RootModel::_root_biomass);
        Internal(LAST_VALUE, &RootModel::_last_value);


        //    external variables
        External(LEAF_DEMAND_SUM, &RootModel::_leaf_demand_sum);
        External(LEAF_LAST_DEMAND_SUM, &RootModel::_leaf_last_demand_sum);
        External(INTERNODE_DEMAND_SUM, &RootModel::_internode_demand_sum);
        External(INTERNODE_LAST_DEMAND_SUM, &RootModel::_internode_last_demand_sum);
        External(PLANT_PHASE, &RootModel::_plant_phase);
        External(PLANT_STATE, &RootModel::_plant_state);
        External(CULM_SURPLUS_SUM, &RootModel::_culm_surplus_sum);
    }

    virtual ~RootModel()
    {}


    void compute(double t, bool /* update */) {
        // parameters
        _P = _parameters.get(t).P;

        // Root Demand Coef
        _root_demand_coef = _coeff1_R_d * std::exp(_coeff2_R_d * (t - _parameters.beginDate + 1)) * (_P * _resp_R_d + 1);

        // Root Demand
        _last_root_demand = _root_demand;
        if (t == _parameters.beginDate) {
            _root_demand = (_leaf_demand_sum + _leaf_last_demand_sum +
                            _internode_demand_sum) * _root_demand_coef;
            _last_value = _root_demand;
            _root_biomass = _root_demand;
        } else {
            if (_plant_state & plant::NOGROWTH) {
                _root_demand = 0;
                _last_value = 0;
            } else {
                if (_plant_phase == plant::ELONG or _plant_phase == plant::PI) {
                    if (_leaf_demand_sum + _leaf_last_demand_sum + _internode_demand_sum
                            + _internode_last_demand_sum == 0) {
                        _root_demand = _last_value * _root_demand_coef;
                    } else {
                        _root_demand = (_leaf_demand_sum + _leaf_last_demand_sum + _internode_demand_sum +
                                        _internode_last_demand_sum) * _root_demand_coef;
                    }
                    if (_leaf_demand_sum + _leaf_last_demand_sum + _internode_demand_sum +
                            _internode_last_demand_sum != 0) {
                        _last_value = _leaf_demand_sum + _leaf_last_demand_sum
                                + _internode_demand_sum + _internode_last_demand_sum;
                    } else {
                        _last_value = 0;
                    }
                    if(_plant_phase == plant::PRE_FLO or _plant_phase == plant::FLO or _plant_phase == plant::END_FILLING) {
                        _root_demand = std::min(_culm_surplus_sum, _root_demand);
                        _surplus = _culm_surplus_sum - _root_demand;
                    } else {
                        _root_demand = _root_demand;
                        _surplus = 0;
                    }
                } else {
                    if (_leaf_demand_sum + _leaf_last_demand_sum == 0) {
                        _root_demand = _last_value * _root_demand_coef;
                    } else {
                        _root_demand = (_leaf_demand_sum + _leaf_last_demand_sum) * _root_demand_coef;
                    }
                    if (_leaf_demand_sum + _internode_demand_sum != 0) {
                        _last_value = _leaf_demand_sum + _internode_demand_sum;
                    } else {
                        _last_value = 0;
                    }
                    _surplus = 0;
                }
            }
            _root_biomass = _root_biomass + _root_demand;
        }
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;

        //    parameters variables
        _coeff1_R_d = _parameters.get("coeff1_R_d");
        _coeff2_R_d = _parameters.get("coeff2_R_d");
        _resp_R_d = _parameters.get("resp_R_d");

        //    computed variables (internal)
        _root_demand = 0;
        _root_demand_coef = 0;
        _last_value = 0;
        _root_biomass = 0;
        _surplus = 0;
        _last_root_demand = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;
    //    parameters
    double _coeff1_R_d;
    double _coeff2_R_d;
    double _resp_R_d;

    //    internals - computed
    double _root_demand_coef;
    double _root_demand;
    double _root_biomass;
    double _last_value;
    double _surplus;
    double _last_root_demand;

    //    externals
    double _P;
    double _leaf_demand_sum;
    double _leaf_last_demand_sum;
    double _internode_demand_sum;
    double _internode_last_demand_sum;
    plant::plant_phase _plant_phase;
    plant::plant_state _plant_state;
    double _culm_surplus_sum;
};

} // namespace model

#endif //ROOT_MODEL_HPP
