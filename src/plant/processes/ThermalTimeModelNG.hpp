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
    enum internals { CULM_BOOL_CROSSED_PLASTO, CULM_BOOL_CROSSED_PHYLLO, CULM_BOOL_CROSSED_LIGULO,
                     CULM_PLASTO_VISU, CULM_PHYLLO_VISU, CULM_LIGULO_VISU,
                     PHENO_STAGE, APP_STAGE, LIG_STAGE, CULM_DD, CULM_DD_PHYLLO, CULM_DD_LIGULO,
                     EDD, EDD_PHYLLO, EDD_LIGULO, TEMP_DD, TEMP_DD_PHYLLO, TEMP_DD_LIGULO };

    enum externals { PLASTO_DELAY, PLASTO, CULM_STOCK, CULM_DEFICIT, DELTA_T,
                     DD, PLASTO_VISU, LIGULO_VISU,
                     IS_FIRST_DAY_OF_INDIVIDUALIZATION, PHYLLO, LIGULO, DD_PHYLLO,
                     DD_LIGULO, PHYLLO_VISU };


    ThermalTimeModelNG() {
        //    computed variables
        Internal(CULM_BOOL_CROSSED_PLASTO, &ThermalTimeModelNG::_culm_bool_crossed_plasto);
        Internal(CULM_BOOL_CROSSED_PHYLLO, &ThermalTimeModelNG::_culm_bool_crossed_phyllo);
        Internal(CULM_BOOL_CROSSED_LIGULO, &ThermalTimeModelNG::_culm_bool_crossed_ligulo);
        Internal(CULM_PLASTO_VISU, &ThermalTimeModelNG::_culm_plastoVisu);
        Internal(CULM_PHYLLO_VISU, &ThermalTimeModelNG::_culm_phylloVisu);
        Internal(CULM_LIGULO_VISU, &ThermalTimeModelNG::_culm_liguloVisu);
        Internal(PHENO_STAGE, &ThermalTimeModelNG::_phenoStage);
        Internal(APP_STAGE, &ThermalTimeModelNG::_appStage);
        Internal(LIG_STAGE, &ThermalTimeModelNG::_ligStage);
        Internal(CULM_DD, &ThermalTimeModelNG::_culm_DD);
        Internal(CULM_DD_PHYLLO, &ThermalTimeModelNG::_culm_DD_phyllo);
        Internal(CULM_DD_LIGULO, &ThermalTimeModelNG::_culm_DD_ligulo);
        Internal(EDD, &ThermalTimeModelNG::_EDD);
        Internal(EDD_PHYLLO, &ThermalTimeModelNG::_EDD_phyllo);
        Internal(EDD_LIGULO, &ThermalTimeModelNG::_EDD_ligulo);
        Internal(TEMP_DD, &ThermalTimeModelNG::_tempDD);
        Internal(TEMP_DD_PHYLLO, &ThermalTimeModelNG::_tempDD_phyllo);
        Internal(TEMP_DD_LIGULO, &ThermalTimeModelNG::_tempDD_ligulo);

        //    external variables
        External(PLASTO_DELAY, &ThermalTimeModelNG::_plasto_delay);
        External(PLASTO, &ThermalTimeModelNG::_plasto);
        External(CULM_STOCK, &ThermalTimeModelNG::_stock);
        External(CULM_DEFICIT, &ThermalTimeModelNG::_deficit);
        External(DELTA_T, &ThermalTimeModelNG::_deltaT);
        External(DD, &ThermalTimeModelNG::_DD);
        External(PLASTO_VISU, &ThermalTimeModelNG::_plastoVisu);
        External(LIGULO_VISU, &ThermalTimeModelNG::_liguloVisu);
        External(IS_FIRST_DAY_OF_INDIVIDUALIZATION, &ThermalTimeModelNG::_is_first_day_of_individualization);
        External(PHYLLO, &ThermalTimeModelNG::_phyllo);
        External(LIGULO, &ThermalTimeModelNG::_ligulo);
        External(DD_PHYLLO, &ThermalTimeModelNG::_DD_phyllo);
        External(DD_LIGULO, &ThermalTimeModelNG::_DD_ligulo);
        External(PHYLLO_VISU, &ThermalTimeModelNG::_phylloVisu);
    }

    virtual ~ThermalTimeModelNG()
    {}

    void compute(double t, bool /* update */) {
        if(_is_first_day_of_individualization) {
            if (_stock + _deficit > 0) {
                _tempDD = _DD + _deltaT + _plasto_delay;
                _tempDD_phyllo = _DD_phyllo + _deltaT + _plasto_delay;
                _tempDD_ligulo = _DD_ligulo + _deltaT + _plasto_delay;

                _culm_bool_crossed_plasto = _tempDD - _plasto;
                _culm_bool_crossed_phyllo = _tempDD_phyllo - _phyllo;
                _culm_bool_crossed_ligulo = _tempDD_ligulo - _ligulo;

                _culm_plastoVisu = _plastoVisu - _plasto_delay;
                _culm_phylloVisu = _phylloVisu - _plasto_delay;
                _culm_liguloVisu = _liguloVisu - _plasto_delay;

                if (_culm_bool_crossed_plasto >= 0) {
                    _EDD = _plasto - _DD;
                    _culm_DD = _tempDD - _plasto;
                    _phenoStage = _phenoStage + 1;
                } else {
                    _EDD = _deltaT + _plasto_delay;
                    _culm_DD = _tempDD;
                }

                if(_culm_bool_crossed_phyllo >= 0) {
                    _EDD_phyllo = _phyllo - _DD_phyllo;
                    _culm_DD_phyllo = _tempDD_phyllo - _phyllo;
                    _appStage = _appStage + 1;
                } else {
                    _EDD_phyllo = _deltaT + _plasto_delay;
                    _culm_DD_phyllo = _tempDD_phyllo;
                }

                if(_culm_bool_crossed_ligulo >= 0) {
                    _EDD_ligulo = _plasto - _DD_ligulo;
                    _culm_DD_ligulo = _tempDD_ligulo - _ligulo;
                    _ligStage = _ligStage + 1;
                } else {
                    _EDD_ligulo = _deltaT + _plasto_delay;
                    _culm_DD_ligulo = _tempDD_ligulo;
                }

            } else {
                _culm_plastoVisu = _plastoVisu + _EDD;
                _culm_phylloVisu = _phylloVisu + _EDD_phyllo;
                _culm_liguloVisu = _liguloVisu + _EDD_ligulo;
            }
        } else {
            if (_stock + _deficit > 0) {
                _tempDD = _culm_DD + _deltaT + _plasto_delay;
                _tempDD_phyllo = _culm_DD_phyllo + _deltaT + _plasto_delay;
                _tempDD_ligulo = _culm_DD_ligulo + _deltaT + _plasto_delay;
                _culm_bool_crossed_plasto = _tempDD - _plasto;
                _culm_bool_crossed_phyllo = _tempDD_phyllo - _phyllo;
                _culm_bool_crossed_ligulo = _tempDD_ligulo - _ligulo;
                _culm_plastoVisu = _culm_plastoVisu - _plasto_delay;
                _culm_phylloVisu = _culm_phylloVisu - _plasto_delay;
                _culm_liguloVisu = _culm_liguloVisu - _plasto_delay;

                if (_culm_bool_crossed_plasto >= 0) {
                    _EDD = _plasto - _culm_DD;
                    _culm_DD = _tempDD - _plasto;
                    _phenoStage = _phenoStage + 1;
                } else {
                    _EDD = _deltaT + _plasto_delay;
                    _culm_DD = _tempDD;
                }

                if (_culm_bool_crossed_phyllo >= 0) {
                    _EDD_phyllo = _phyllo - _culm_DD_phyllo;
                    _culm_DD_phyllo = _tempDD_phyllo - _phyllo;
                    _appStage = _appStage + 1;
                } else {
                    _EDD_phyllo = _deltaT + _plasto_delay;
                    _culm_DD_phyllo = _tempDD_phyllo;
                }

                if (_culm_bool_crossed_ligulo >= 0) {
                    _EDD_ligulo = _ligulo - _culm_DD_ligulo;
                    _culm_DD_ligulo = _tempDD_ligulo - _ligulo;
                    _ligStage = _ligStage + 1;
                } else {
                    _EDD_ligulo = _deltaT + _plasto_delay;
                    _culm_DD_ligulo = _tempDD_ligulo;
                }
            } else {
                _culm_plastoVisu = _culm_plastoVisu + _EDD;
                _culm_phylloVisu = _phylloVisu + _EDD_phyllo;
                _culm_liguloVisu = _culm_liguloVisu + _EDD_ligulo;
            }
        }
    }



    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        //    paramaters variables
        _plasto_init = _parameters.get("plasto_init");
        _phyllo_init = _parameters.get("phyllo_init");
        _ligulo_init = _parameters.get("ligulo_init");

        //    computed variables
        _culm_bool_crossed_plasto = 0;
        _culm_bool_crossed_phyllo = 0;
        _culm_bool_crossed_ligulo = 0;
        _culm_plastoVisu = _plasto_init;
        _culm_phylloVisu = _phyllo_init;
        _culm_liguloVisu = _ligulo_init;
        _phenoStage = 0;
        _appStage = 0;
        _ligStage = 0;
        _culm_DD = 0;
        _EDD = 0;
        _tempDD = 0;
        _culm_DD_phyllo = 0;
        _EDD_phyllo = 0;
        _tempDD_phyllo = 0;
        _culm_DD_ligulo = 0;
        _EDD_ligulo = 0;
        _tempDD_ligulo = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;
    //    parameters
    double _plasto_init;
    double _phyllo_init;
    double _ligulo_init;

    //    internals
    double _culm_bool_crossed_plasto;
    double _culm_bool_crossed_phyllo;
    double _culm_bool_crossed_ligulo;

    double _culm_plastoVisu;
    double _culm_phylloVisu;
    double _culm_liguloVisu;

    int _phenoStage;
    double _culm_DD;
    double _EDD;
    double _tempDD;

    int _appStage;
    double _culm_DD_phyllo;
    double _EDD_phyllo;
    double _tempDD_phyllo;

    int _ligStage;
    double _culm_DD_ligulo;
    double _EDD_ligulo;
    double _tempDD_ligulo;


    //    externals
    double _plasto;
    double _phyllo;
    double _ligulo;

    double _plasto_delay;

    double _stock;
    double _deficit;
    double _deltaT;

    double _DD;
    double _DD_phyllo;
    double _DD_ligulo;

    double _plastoVisu;
    double _phylloVisu;
    double _liguloVisu;

    bool _is_first_day_of_individualization;

};

} // namespace model
#endif
