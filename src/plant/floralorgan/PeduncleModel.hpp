/**
 * @file ecomeristem/peduncle/Model.hpp
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

#include <defines.hpp>

using namespace std;

namespace model {

class PeduncleModel : public AtomicModel < PeduncleModel >
{
public:

    enum internals { IS_MATURE, LENGTH_PREDIM, DIAMETER_PREDIM, REDUCTION_INER, INER,
                     LENGTH, VOLUME, EXP_TIME, BIOMASS, DEMAND, LAST_DEMAND, FIRST_DAY };
    enum externals { PLANT_PHASE, CULM_PHASE, INTER_PREDIM, INTER_DIAM, FTSW,
                     EDD, DELTA_T, PLASTO, LIGULO, FCSTR, FCSTRI, TEST_IC };


    PeduncleModel(int index, bool is_on_mainstem):
        _index(index),
        _is_on_mainstem(is_on_mainstem)
    {
        Internal(IS_MATURE, &PeduncleModel::_is_mature);
        Internal(LENGTH_PREDIM, &PeduncleModel::_length_predim);
        Internal(DIAMETER_PREDIM, &PeduncleModel::_diameter_predim);
        Internal(REDUCTION_INER, &PeduncleModel::_reduction_iner);
        Internal(INER, &PeduncleModel::_iner);
        Internal(LENGTH, &PeduncleModel::_length);
        Internal(VOLUME, &PeduncleModel::_volume);
        Internal(EXP_TIME, &PeduncleModel::_exp_time);
        Internal(BIOMASS, &PeduncleModel::_biomass);
        Internal(DEMAND, &PeduncleModel::_demand);
        Internal(LAST_DEMAND, &PeduncleModel::_last_demand);
        Internal(FIRST_DAY, &PeduncleModel::_first_day);

        External(PLANT_PHASE, &PeduncleModel::_plant_phase);
        External(CULM_PHASE, &PeduncleModel::_culm_phase);
        External(INTER_PREDIM, &PeduncleModel::_inter_predim); //FirstNonVegetativeInternode
        External(INTER_DIAM, &PeduncleModel::_inter_diam); //FirstNonVegetativeInternode
        External(FTSW, &PeduncleModel::_ftsw);
        External(EDD, &PeduncleModel::_edd);
        External(DELTA_T, &PeduncleModel::_delta_t);
        External(PLASTO, &PeduncleModel::_plasto);
        External(LIGULO, &PeduncleModel::_ligulo);
        External(FCSTR, &PeduncleModel::_fcstr);
        External(FCSTRI, &PeduncleModel::_fcstrI);
        External(TEST_IC, &PeduncleModel::_test_ic);
    }

    virtual ~PeduncleModel()
    { }

    void compute(double t, bool /* update */)
    {
        _p = _parameters.get(t).P;
        if(t == _first_day) {
            //Peduncle Length Predim
            _length_predim = _ratio_in_ped * _inter_predim;

            //Peduncle Diameter Predim
            _diameter_predim = _peduncle_diam * _inter_diam;
        }

        if(!_is_mature and _plant_phase != plant::MATURITY and _culm_phase != culm::FLO) {

            //ReductionINER
            if(_wbmodel == 2) {
                _reduction_iner = std::max(1e-4, (std::min(1.,_fcstrI * (1. + (_p * _respINER)))) * _test_ic);
            } else {
                if (_ftsw < _thresINER) {
                    _reduction_iner = std::max(1e-4, ((1./_thresINER) * _ftsw) * (1. + (_p * _respINER))) * _test_ic;
                } else {
                    _reduction_iner = 1. + _p * _respINER * _test_ic;
                }
            }

            _iner = _length_predim * _reduction_iner / (_phenostage_pre_flo_to_flo * _ligulo);

            //Length and Exp time
            if (t == _first_day) {
                _length = _iner * _edd;
                _exp_time = (_length_predim - _length) / _iner;
            } else {
                _exp_time = (_length_predim - _length) / _iner;
                _length = std::min(_length_predim, _length + (_iner * std::min(_delta_t, _exp_time)));
            }

            //Volume
            double radius = _diameter_predim / 2;
            _volume = _length * 3.141592653589793238462643383280 * radius * radius;

            //Biomass and Demand
            if(t == _first_day) {
                _biomass = _volume * _density;
                _demand = _biomass;
            } else {
                _demand = (_density * _volume) - _biomass;
                _biomass = _density * _volume;
            }

            if(_length == _length_predim) {
                _is_mature = true;
                _last_demand = _demand;
                _demand = 0;
            }
        }
        _demand = 0;
        _last_demand = _demand;
    }

    void init(double t, const ecomeristem::ModelParameters&  parameters )
    {
        _parameters = parameters;

        // parameters
        _ratio_in_ped = _parameters.get("ratio_INPed");
        _peduncle_diam = _parameters.get("peduncle_diam");
        _thresINER = _parameters.get("thresINER");
        _respINER = _parameters.get("resp_LER");
        _density = _parameters.get("density_IN2");
        _wbmodel = _parameters.get("wbmodel");
        _phenostage_pre_flo_to_flo = parameters.get("phenostage_PRE_FLO_to_FLO");

        // internals
        _is_mature = false;
        _length_predim = 0;
        _diameter_predim = 0;
        _reduction_iner = 0;
        _iner = 0;
        _length = 0;
        _volume = 0;
        _exp_time = 0;
        _biomass = 0;
        _first_day = t;
        _demand = 0;
        _last_demand = 0;
    }
private:
    ecomeristem::ModelParameters  _parameters;

    // attributes
    int _index;
    bool _is_first_internode;
    bool _is_on_mainstem;

    // parameters
    double _ratio_in_ped;
    double _peduncle_diam;
    double _p;
    double _thresINER;
    double _respINER;
    double _density;
    double _wbmodel;
    double _phenostage_pre_flo_to_flo;

    // internals
    bool _is_mature;
    double _length_predim;
    double _diameter_predim;
    double _reduction_iner;
    double _iner;
    double _length;
    double _volume;
    double _exp_time;
    double _biomass;
    double _demand;
    double _last_demand;
    double _first_day;

    // externals
    plant::plant_phase _plant_phase;
    culm::culm_phase _culm_phase;
    double _inter_predim;
    double _inter_diam;
    double _ftsw;
    double _edd;
    double _delta_t;
    double _plasto;
    double _ligulo;
    double _fcstr;
    double _fcstrI;
    double _test_ic;

};

}
