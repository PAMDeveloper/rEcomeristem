/**
 * @file ecomeristem/panicle/Model.hpp
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

#define _USE_MATH_DEFINES
#include <math.h>

namespace model {

class PanicleModel : public AtomicModel < PanicleModel >
{
public:
    enum internals { PREFLO_PASSED, GRAIN_NB, FERTILE_GRAIN_NB, RESERVOIR_DISPO,
                     DAY_DEMAND, WEIGHT, FILLED_GRAIN_NB };
    enum externals { DELTA_T, FCSTR, TEST_IC, PLANT_PHASE };

    PanicleModel()
    {
        Internal(PREFLO_PASSED, &PanicleModel::_preflo_passed);
        Internal(GRAIN_NB, &PanicleModel::_grain_nb);
        Internal(FERTILE_GRAIN_NB, &PanicleModel::_fertile_grain_nb);
        Internal(RESERVOIR_DISPO, &PanicleModel::_reservoir_dispo);
        Internal(DAY_DEMAND, &PanicleModel::_day_demand);
        Internal(WEIGHT, &PanicleModel::_weight);
        Internal(FILLED_GRAIN_NB, &PanicleModel::_filled_grain_nb);

        External(DELTA_T, &PanicleModel::_delta_t);
        External(FCSTR, &PanicleModel::_fcstr);
        External(TEST_IC, &PanicleModel::_test_ic);
        External(PLANT_PHASE, &PanicleModel::_plant_phase);
    }

    virtual ~PanicleModel()
    { }

    void compute(double t, bool /* update */)
    {
        // Floraison
        if (_plant_phase == plant::FLO) {
            if (!_preflo_passed) {
                // Panicle Fertile Grain Nb
                _fertile_grain_nb = _fertile_grain_nb * _fcstr * _test_ic;
                _preflo_passed = true;
            }

            // Panicle Reservoir Dispo
            _reservoir_dispo = std::max(0., _fertile_grain_nb * _gdw - _weight);

            // Panicle Day Demand
            _day_demand = std::min(_grain_filling_rate * _delta_t * _fcstr * _test_ic, _reservoir_dispo);

            // Panicle Weight
            _weight = _weight + _day_demand;

            // Panicle Filled Grain Nb
            _filled_grain_nb = (_fertile_grain_nb * _weight) / (_fertile_grain_nb * _gdw);
        } else {
            _day_demand = 0;
        }

        if(!_preflo_passed) {
            _grain_nb = _grain_nb + (_spike_creation_rate * _delta_t * _fcstr * _test_ic);
            _fertile_grain_nb = _grain_nb;
        }

    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        // Parameters
        _spike_creation_rate = _parameters.get("spike_creation_rate");
        _grain_filling_rate = _parameters.get("grain_filling_rate");
        _gdw = _parameters.get("gdw");

        // Internals
        _preflo_passed = false;
        _grain_nb = 0;
        _fertile_grain_nb = 0;;
        _reservoir_dispo = 0;
        _day_demand = 0;
        _weight = 0;
        _filled_grain_nb = 0;
    }


private:
    ecomeristem::ModelParameters _parameters;

    // Parameters
    double _spike_creation_rate;
    double _grain_filling_rate;
    double _gdw;

    // Internals
    bool _preflo_passed;
    double _grain_nb;
    double _fertile_grain_nb;
    double _reservoir_dispo;
    double _day_demand;
    double _weight;
    double _filled_grain_nb;

    // Externals
    double _delta_t;
    double _fcstr;
    double _test_ic;
    plant::plant_phase _plant_phase;

};

} //namespace
