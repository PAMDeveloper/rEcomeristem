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
        // _etp = computeETP();
        _water_supply = _parameters.get(t).Irrigation;
        if(wbmodel == 1) {
            //Pot Waterbalance model [phenoarch 2017]

            //FTSW
            _ftsw = (_swc-pf) / RU1;

            //cstr
            _cstr = (_ftsw < ThresTransp) ? std::max(1e-4, _ftsw * 1. / ThresTransp) : 1;

            //fcstr
            _fcstr = std::sqrt(_cstr);

            //See FAO ETP http://www.fao.org/3/X0490E/x0490e0c.htm#chapter%207
            if(_De < REW) { //REW and TEW are soil parameters
                _Kr = 1;
            } else {
                _Kr = (TEW - _De)/(TEW - REW);
            }
            _Ke = _Kr * (Ke_init); //with _Ke_init = Kcmax - Kcpot; estimate on data

            //transpiration
            _transpiration = std::min(_swc, (Kcpot * std::min(_etp, ETPmax) * _interc * _cstr) / Density);

            //evaporation
            _evaporation = (_Ke * _etp) * 0.0283;

            //SWC
            _swc = _swc - _transpiration - _evaporation + _water_supply;

            //De
            _De = _De + _evaporation - _water_supply;


        } else {
            _water_supply = _parameters.get(t).Irrigation;

            //Field waterbalance model [BFF 2014-2015-2016]
            _cstr = 1;
            _transpiration = _swc;
            _ftsw = 1;
            if(_water_supply == 1) {
                _stressdays = std::min(10.,_stressdays + 1);
                _psib = (pot/10)*_stressdays;
                _fcstrA = std::min(1.,std::max(0.,(((thresAssim+thresLER)+(stressBP2+stressBP))-_psib)/(((thresAssim+thresLER)+(stressBP2+stressBP))-(stressBP2+stressBP))));
                _fcstrL = std::min(1.,std::max(0.,((thresLER+stressBP)-_psib)/((thresLER+stressBP)-stressBP)));
                _fcstrI = std::min(1.,std::max(0.,((thresINER+stressBP)-_psib)/((thresINER+stressBP)-stressBP)));
                _fcstrLlen = std::min(1.,std::max(0.,((thresLEN+stressBP)-_psib)/((thresLEN+stressBP)-stressBP)));
            } else {
                _fcstr = 1;
                _fcstrA = 1;
                _fcstrI = 1;
                _fcstrL = 1;
                _stressdays = 0;
            }
        }
    }

    //double computeKe() {
       // TEW = 1000*(1.7-(0.5*0.34)*0.03);
       // if(true) {
       //     Kr = 1;
       // } else {
       //     Kr = (TEW - De)/(TEW - REW);
      //  }
      //  Ke = Kr(Kcmax - Kcb)
    //}

    double computeETP(double t) {
        // Weather and param variables TODO
        double RgMax;
        double RgCalc;
        double TMin;
        double TMax;
        double HMin;
        double HMax;
        double HMoyCalc;
        double TMoyCalc;
        double Vt;
        double Altitude;
        double TMoyPrec;
        double VPDCalc;

        double eActual; double eSat; double RgRgMax; double TLat; double delta; double KPsy; double Eaero; double Erad; double Rn; double G; double ETo;
        eSat = 0.3054 * (exp(17.27 * TMax * 1.0 / (TMax + 237.3)) + exp(17.27 * TMin * 1.0 / (TMin + 237.3)));
        if ((HMax == 0)) {
            eActual = eSat * HMoyCalc * 1.0 / 100;
        }
        else {
            eActual = 0.3054 * (exp(17.27 * TMax * 1.0 / (TMax + 237.3)) * HMin * 1.0 / 100 + exp(17.27 * TMin * 1.0 / (TMin + 237.3)) * HMax * 1.0 / 100);
        }

        // VPD calc or weather file ?
        VPDCalc = eSat - eActual;
        RgRgMax = min(1.0,RgCalc * 1.0 / RgMax);
        Rn = 0.77 * RgCalc - (1.35 * RgRgMax - 0.35) * (0.34 - 0.14 * std::pow(eActual, 0.5)) * (pow(TMax + 273.16, 4) + std::pow(TMin + 273.16, 4)) * 2.45015 * std::pow(10, -9);

        // chaleur latente de vaporisation de l'eau
        TLat = 2.501 - 2.361 * std::pow(10, -3) * TMoyCalc;
        //  pente de la courbe de pression de vapeur saturante en kPa/°C
        delta = 4098 * (0.6108 * exp(17.27 * TMoyCalc * 1.0 / (TMoyCalc + 237.3))) * 1.0 / std::pow(TMoyCalc + 237.3, 2);
        // constante psychrométrique en kPa/°C
        KPsy = 0.00163 * 101.3 * std::pow(1 - (0.0065 * Altitude * 1.0 / 293), 5.26) * 1.0 / TLat;
        // Radiative
        G = 0.38 * (TMoyCalc - TMoyPrec);
        Erad = 0.408 * (Rn - G) * delta * 1.0 / (delta + KPsy * (1 + 0.34 * Vt));
        // Partie évaporative de ET0 = Eaéro
        Eaero = (900 * 1.0 / (TMoyCalc + 273.16)) * ((eSat - eActual) * Vt) * KPsy * 1.0 / (delta + KPsy * (1 + 0.34 * Vt));
        ETo = Erad + Eaero;
        return(ETo);
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
        stressBP2 = parameters.get("stressBP2");
        pot = parameters.get("psib");
        pf = parameters.get("pf");;
        swc_init = parameters.get("swc_init");
        wbmodel = parameters.get("wbmodel");
        stressD = parameters.get("stressD");
        TEW = parameters.get("TEW");
        REW = parameters.get("REW");
        Ke_init = parameters.get("Ke_init");

        //    computed variables
        _cstr = 1;
        _fcstr = 1;
        _ftsw = 1;
        _swc = swc_init;
        _transpiration = 0;
        _psib = 0;
        _stressdays = 0;
        _fcstrA = 1;
        _fcstrL = 1;
        _fcstrI = 1;
        _fcstrLlen = 1;
        _De = 0;
        _evaporation = 0;
        _Ke = 0;
        _Kr = 1;

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
    double wbmodel;
    double stressBP;
    double stressBP2;
    //Pot WB
    double pf;
    double swc_init;
    double stressD;
    double TEW;
    double REW;
    double Ke_init;

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
    //Pot WB
    double _De;
    double _evaporation;
    double _Ke;
    double _Kr;

    //  externals
    double _interc;
};

} // namespace model
#endif
