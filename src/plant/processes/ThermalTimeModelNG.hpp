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


#ifndef THERMAL_TIME_MODEL_NG_HPP
#define THERMAL_TIME_MODEL_NG_HPP

#include <defines.hpp>

namespace model {

class ThermalTimeModelNG : public AtomicModel < ThermalTimeModelNG >
{
public:
    enum internals { CULM_BOOL_CROSSED_PLASTO, CULM_PLASTO_VISU, CULM_LIGULO_VISU,
                     PHENO_STAGE, CULM_DD, EDD, TEMP_DD };

    enum externals { PLASTO_DELAY, PLASTO, CULM_STOCK, CULM_DEFICIT, DELTA_T,
                     BOOL_CROSSED_PLASTO, DD, PLASTO_VISU, LIGULO_VISU,
                     IS_FIRST_DAY_OF_INDIVIDUALIZATION };


    ThermalTimeModelNG() {
        //    computed variables
        Internal(CULM_BOOL_CROSSED_PLASTO, &ThermalTimeModelNG::_culm_bool_crossed_plasto);
        Internal(CULM_PLASTO_VISU, &ThermalTimeModelNG::_culm_plastoVisu);
        Internal(CULM_LIGULO_VISU, &ThermalTimeModelNG::_culm_liguloVisu);
        Internal(PHENO_STAGE, &ThermalTimeModelNG::_phenoStage);
        Internal(CULM_DD, &ThermalTimeModelNG::_culm_DD);
        Internal(EDD, &ThermalTimeModelNG::_EDD);
        Internal(TEMP_DD, &ThermalTimeModelNG::_tempDD);

        //    external variables
        External(PLASTO_DELAY, &ThermalTimeModelNG::_plasto_delay);
        External(PLASTO, &ThermalTimeModelNG::_plasto);
        External(CULM_STOCK, &ThermalTimeModelNG::_stock);
        External(CULM_DEFICIT, &ThermalTimeModelNG::_deficit);
        External(DELTA_T, &ThermalTimeModelNG::_deltaT);
        External(BOOL_CROSSED_PLASTO, &ThermalTimeModelNG::_bool_crossed_plasto);
        External(DD, &ThermalTimeModelNG::_DD);
        External(PLASTO_VISU, &ThermalTimeModelNG::_plastoVisu);
        External(LIGULO_VISU, &ThermalTimeModelNG::_liguloVisu);
        External(IS_FIRST_DAY_OF_INDIVIDUALIZATION, &ThermalTimeModelNG::_is_first_day_of_individualization);

    }

    virtual ~ThermalTimeModelNG()
    {}



    void compute(double t, bool /* update */) {
        if(_is_first_day_of_individualization) {
            if (_stock + _deficit > 0) {
                _tempDD = _DD + _deltaT + _plasto_delay;
                _culm_bool_crossed_plasto = _tempDD - _plasto;
                _culm_plastoVisu = _plastoVisu - _plasto_delay;
                _culm_liguloVisu = _liguloVisu - _plasto_delay;

                if (_culm_bool_crossed_plasto >= 0) {
                    _EDD = _plasto - _DD;
                    _culm_DD = _tempDD - _plasto;
                    _phenoStage = _phenoStage + 1;
                } else {
                    _EDD = _deltaT + _plasto_delay;
                    _culm_DD = _tempDD;
                }
            } else {
                _culm_plastoVisu = _plastoVisu + _EDD;
                _culm_liguloVisu = _liguloVisu + _EDD;
            }
        } else {
            if (_stock + _deficit > 0) {
                _tempDD = _culm_DD + _deltaT + _plasto_delay;

                _culm_bool_crossed_plasto = _tempDD - _plasto;
                _culm_plastoVisu = _culm_plastoVisu - _plasto_delay;
                _culm_liguloVisu = _culm_liguloVisu - _plasto_delay;

                if (_culm_bool_crossed_plasto >= 0) {
                    _EDD = _plasto - _culm_DD;
                    _culm_DD = _tempDD - _plasto;
                    _phenoStage = _phenoStage + 1;
                } else {
                    _EDD = _deltaT + _plasto_delay;
                    _culm_DD = _tempDD;
                }
            } else {
                _culm_plastoVisu = _culm_plastoVisu + _EDD;
                _culm_liguloVisu = _culm_liguloVisu + _EDD;
            }
        }
    }



    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        //    paramaters variables
        _coef_ligulo = _parameters.get("coef_ligulo1");
        _plasto_init = _parameters.get("plasto_init");

        //    computed variables
        _culm_bool_crossed_plasto = 0;
        _culm_plastoVisu = _plasto_init;
        _culm_liguloVisu = _plasto_init * _coef_ligulo;
        _phenoStage = 0;
        _culm_DD = 0;
        _EDD = 0;
        _tempDD = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;
    //    parameters
    double _coef_ligulo;
    double _plasto_init;

    //    internals
    double _culm_bool_crossed_plasto;
    double _culm_plastoVisu;
    double _culm_liguloVisu;
    int _phenoStage;
    double _culm_DD;
    double _EDD;
    double _tempDD;

    //    externals
    double _plasto;
    double _plasto_delay;
    double _stock;
    double _deficit;
    double _deltaT;
    double _bool_crossed_plasto;
    double _DD;
    double _plastoVisu;
    double _liguloVisu;
    bool _is_first_day_of_individualization;

};

} // namespace model
#endif
