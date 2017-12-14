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


#ifndef THERMAL_TIME_MODEL_HPP
#define THERMAL_TIME_MODEL_HPP

#include <defines.hpp>

namespace model {

class ThermalTimeModel : public AtomicModel < ThermalTimeModel >
{
public:
    enum internals { BOOL_CROSSED_PLASTO, BOOL_CROSSED_PHYLLO, BOOL_CROSSED_LIGULO,
                     PLASTO_VISU, PHYLLO_VISU, LIGULO_VISU, PHENO_STAGE, APP_STAGE, LIG_STAGE,
                     SLA, TEMP_DD, TEMP_DD_PHYLLO, TEMP_DD_LIGULO, DD, EDD, DD_PHYLLO, DD_LIGULO, EDD_PHYLLO, EDD_LIGULO };

    enum externals {  DELTA_T, PLASTO_DELAY, PLASTO, PHYLLO, LIGULO, STOCK };


    ThermalTimeModel() {
        //    computed variables
        Internal(BOOL_CROSSED_PLASTO, &ThermalTimeModel::_boolCrossedPlasto);
        Internal(BOOL_CROSSED_PHYLLO, &ThermalTimeModel::_boolCrossedPhyllo);
        Internal(BOOL_CROSSED_LIGULO, &ThermalTimeModel::_boolCrossedLigulo);
        Internal(PLASTO_VISU, &ThermalTimeModel::_plastoVisu);
        Internal(PHYLLO_VISU, &ThermalTimeModel::_phylloVisu);
        Internal(LIGULO_VISU, &ThermalTimeModel::_liguloVisu);
        Internal(PHENO_STAGE, &ThermalTimeModel::_phenoStage);
        Internal(APP_STAGE, &ThermalTimeModel::_appStage);
        Internal(LIG_STAGE, &ThermalTimeModel::_ligStage);
        Internal(TEMP_DD, &ThermalTimeModel::_tempDD);
        Internal(DD, &ThermalTimeModel::_DD);
        Internal(EDD, &ThermalTimeModel::_EDD);
        Internal(TEMP_DD_PHYLLO, &ThermalTimeModel::_tempDD_phyllo);
        Internal(DD_PHYLLO, &ThermalTimeModel::_DD_phyllo);
        Internal(EDD_PHYLLO, &ThermalTimeModel::_EDD_phyllo);
        Internal(TEMP_DD_LIGULO, &ThermalTimeModel::_tempDD_ligulo);
        Internal(DD_LIGULO, &ThermalTimeModel::_DD_ligulo);
        Internal(EDD_LIGULO, &ThermalTimeModel::_EDD_ligulo);

        //    external variables
        External(PLASTO_DELAY, &ThermalTimeModel::_plasto_delay);
        External(PLASTO, &ThermalTimeModel::_plasto);
        External(PHYLLO, &ThermalTimeModel::_phyllo);
        External(LIGULO, &ThermalTimeModel::_ligulo);
        External(STOCK, &ThermalTimeModel::_stock);
        External(DELTA_T, &ThermalTimeModel::_deltaT);
    }

    virtual ~ThermalTimeModel()
    {}



    void compute(double t, bool /* update */) {
        if (_stock != 0) {
            _tempDD = _DD + _deltaT + _plasto_delay;
            _tempDD_phyllo = _DD_phyllo + _deltaT + _plasto_delay;
            _tempDD_ligulo = _DD_ligulo + _deltaT + _plasto_delay;

            _boolCrossedPlasto = _tempDD - _plasto;
            _boolCrossedPhyllo = _tempDD_phyllo - _phyllo;
            _boolCrossedLigulo = _tempDD_ligulo - _ligulo;

            _plastoVisu = _plastoVisu - _plasto_delay;
            _phylloVisu = _phylloVisu - _plasto_delay;
            _liguloVisu = _liguloVisu - _plasto_delay;

            if (_boolCrossedPlasto >= 0) {
                _EDD = _plasto - _DD;
                _DD = _tempDD - _plasto;
                _phenoStage = _phenoStage + 1;
            } else {
                _EDD = _deltaT + _plasto_delay;
                _DD = _tempDD;
            }

            if (_boolCrossedPhyllo >= 0) {
                _EDD_phyllo = _phyllo - _DD_phyllo;
                _DD_phyllo = _tempDD_phyllo - _phyllo;
                _appStage = _appStage + 1;
            } else {
                _EDD_phyllo = _deltaT + _plasto_delay;
                _DD_phyllo = _tempDD_phyllo;
            }

            if (_boolCrossedLigulo >= 0) {
                _EDD_ligulo = _ligulo - _DD_ligulo;
                _DD_ligulo = _tempDD_ligulo - _ligulo;
                _ligStage = _ligStage + 1;
            } else {
                _EDD_ligulo = _deltaT + _plasto_delay;
                _DD_ligulo = _tempDD_ligulo;
            }

        } else {
            _plastoVisu = _plastoVisu + _EDD;
            _phylloVisu = _phylloVisu + _EDD;
            _liguloVisu = _liguloVisu + _EDD;
        }
    }



    void init(double /*t*/, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        //    paramaters variables
        _plasto_init = _parameters.get("plasto_init");
        _phyllo_init = _parameters.get("phyllo_init");
        _ligulo_init = _parameters.get("ligulo_init");


        //    computed variables
        _boolCrossedPlasto = 0;
        _boolCrossedPhyllo = 0;
        _boolCrossedLigulo = 0;

        _plastoVisu = _plasto_init;
        _phylloVisu = _phyllo_init;
        _liguloVisu = _ligulo_init;

        _phenoStage = 4;
        _appStage = 1;
        _ligStage = 0;

        _DD = 0;
        _EDD = 0;
        _DD_phyllo = 0;
        _EDD_phyllo = 0;
        _DD_ligulo = 0;
        _EDD_ligulo = 0;

        _tempDD = 0;
        _tempDD_phyllo = 0;
        _tempDD_ligulo = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;
    //    parameters
    double _plasto_init;
    double _phyllo_init;
    double _ligulo_init;

    //    internals
    double _boolCrossedPlasto;
    double _boolCrossedPhyllo;
    double _boolCrossedLigulo;

    double _plastoVisu;
    double _phylloVisu;
    double _liguloVisu;

    int _phenoStage;
    int _appStage;
    int _ligStage;

    double _DD;
    double _EDD;
    double _DD_phyllo;
    double _EDD_phyllo;
    double _DD_ligulo;
    double _EDD_ligulo;

    double _tempDD;
    double _tempDD_phyllo;
    double _tempDD_ligulo;

    //    externals
    double _plasto;
    double _phyllo;
    double _ligulo;
    double _plasto_delay;
    double _stock;
    double _deltaT;
};

} // namespace model
#endif
