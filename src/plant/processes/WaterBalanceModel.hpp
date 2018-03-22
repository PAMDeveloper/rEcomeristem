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


#ifndef WATER_BALANCE_MODEL_HPP
#define WATER_BALANCE_MODEL_HPP

#include <defines.hpp>

namespace model {

class WaterBalanceModel : public AtomicModel < WaterBalanceModel >
{
public:
    enum internals { CSTR, FCSTR, FCSTRA, FCSTRI, FCSTRL, FCSTRLLEN, FTSW, TRANSPIRATION, SWC, PSIB };
    enum externals { INTERC };


    WaterBalanceModel() {
        //    computed variables
        Internal(CSTR, &WaterBalanceModel::_cstr);
        Internal(FCSTR, &WaterBalanceModel::_fcstr);
        Internal(FCSTRA, &WaterBalanceModel::_fcstrA);
        Internal(FCSTRI, &WaterBalanceModel::_fcstrI);
        Internal(FCSTRL, &WaterBalanceModel::_fcstrL);
        Internal(FCSTRLLEN, &WaterBalanceModel::_fcstrLlen);
        Internal(FTSW, &WaterBalanceModel::_ftsw);
        Internal(TRANSPIRATION, &WaterBalanceModel::_transpiration);
        Internal(SWC, &WaterBalanceModel::_swc);
        Internal(PSIB, &WaterBalanceModel::_psib);

        External(INTERC, &WaterBalanceModel::_interc);
    }

    virtual ~WaterBalanceModel()
    {}


    void compute(double t, bool /* update */) {
        //parameters
        _etp = _parameters.get(t).Etp;
        _water_supply = 0;
        if(_wbmodel == 1) {
            //Waterbalance model
            //FTSW
            _ftsw = _swc / RU1;

            //cstr
            _cstr = (_ftsw < ThresTransp) ? std::max(1e-4, _ftsw * 1. / ThresTransp) : 1;

            //fcstr
            _fcstr = std::sqrt(_cstr);

            //transpiration
            _transpiration = std::min(_swc, (Kcpot * std::min(_etp, ETPmax) * _interc * _cstr) / Density);

            //SWC
            _swc = _swc - _transpiration + _water_supply;
        } else {
            _water_supply = _parameters.get(t).Irrigation;

            //Field waterbalance model
            _cstr = 1;
            _transpiration = _swc;
            _ftsw = 1;
            if(_water_supply == 1) {
                _stressdays = std::min(10.,_stressdays + 1);
                _psib = (pot/10)*_stressdays;
                _fcstrA = std::min(1.,std::max(0.,(stressBP-_psib)/(stressBP-thresAssim)));
                _fcstrL = std::min(1.,std::max(0.,(stressBP-_psib)/(stressBP-thresLER)));
                _fcstrI = std::min(1.,std::max(0.,(stressBP-_psib)/(stressBP-thresINER)));
                _fcstrLlen = std::min(1.,std::max(0.,(stressBP-_psib)/(stressBP-thresLEN)));
            } else {
                _fcstr = 1;
                _fcstrA = 1;
                _fcstrI = 1;
                _fcstrL = 1;
                _stressdays = 0;
            }
        }
    }


    void init(double /*t*/, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;

        //    paramaters variables
        ThresTransp = parameters.get("thresTransp");
        RU1 = parameters.get("RU1");
        ETPmax = parameters.get("ETPmax");
        Kcpot = parameters.get("Kcpot");
        Density = parameters.get("density");
        thresLER = parameters.get("thresLER");
        thresINER = parameters.get("thresINER");
        thresAssim = parameters.get("thresAssim");
        thresLEN = parameters.get("thresLEN");
        stressBP = parameters.get("stressBP");
        pot = parameters.get("psib");

        _wbmodel = parameters.get("wbmodel");
        //    computed variables
        _cstr = 1;
        _fcstr = 1;
        _ftsw = 1;
        _swc = RU1;
        _transpiration = 0;
        _psib = 0;
        _stressdays = 0;
        _fcstrA = 1;
        _fcstrL = 1;
        _fcstrI = 1;
        _fcstrLlen = 1;
    }

private:
    ecomeristem::ModelParameters _parameters;
    // parameters
    double ETPmax;
    double Kcpot;
    double Density;
    double RU1;
    double ThresTransp;
    //Field WB
    double thresLER;
    double thresAssim;
    double thresINER;
    double thresLEN;
    double pot;
    double _wbmodel;
    double stressBP;
    //meteo
    double _etp;
    double _water_supply;

    //    internals (computed)
    double _transpiration;
    double _ftsw;
    double _swc;
    double _cstr;
    double _fcstr;
    double _fcstrA;
    double _fcstrL;
    double _fcstrI;
    double _fcstrLlen;
    //Field WB
    double _psib;
    double _stressdays;

    //  externals
    double _interc;
};

} // namespace model
#endif
