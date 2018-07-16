/**
 * @file ecomeristem/plant/stock/Model.hpp
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

// GECROS MODEL, see "Crop Systems Biology" for equations


#ifndef INTERCEPTION_MODEL_HPP
#define INTERCEPTION_MODEL_HPP

#include <defines.hpp>
#include <random>

namespace model {

class InterceptionModel : public AtomicModel < InterceptionModel >
{
public:
    enum internals { INTERC, RG, DECLINATION, AHS, SUNRISE, SUNSET, DAYDURATION, MEAN,
                     SIGMA, SCD, SINLD, COSLD, DAYLENGTH, SOLARHEIGHT, REXT, TRANSATM,
                     PROPRGRD, KPD, LINTERC, LAI, KPD_15, KPD_45, KPD_75, ALPHA, DOY};

    enum externals { PAI };


    InterceptionModel() {
        //  computed variables
        Internal(INTERC, &InterceptionModel::_interc);
        Internal(RG, &InterceptionModel::_Rg);
        Internal(DECLINATION, &InterceptionModel::_declination);
        Internal(AHS, &InterceptionModel::_ahs);
        Internal(SUNRISE, &InterceptionModel::_sunrise);
        Internal(SUNSET, &InterceptionModel::_sunset);
        Internal(DAYDURATION, &InterceptionModel::_dayDuration);
        Internal(MEAN, &InterceptionModel::_mean);
        Internal(SIGMA, &InterceptionModel::_sigma);
        Internal(SCD, &InterceptionModel::_scd);
        Internal(SINLD, &InterceptionModel::_sinLD);
        Internal(COSLD, &InterceptionModel::_cosLD);
        Internal(DAYLENGTH, &InterceptionModel::_dayLength);
        Internal(SOLARHEIGHT, &InterceptionModel::_solarHeight);
        Internal(REXT, &InterceptionModel::_rExt);
        Internal(TRANSATM, &InterceptionModel::_TransAtm);
        Internal(PROPRGRD, &InterceptionModel::_propRgRd);
        Internal(KPD, &InterceptionModel::_kpd);
        Internal(LINTERC, &InterceptionModel::_Linterc);
        Internal(KPD_15, &InterceptionModel::_kpd_15);
        Internal(KPD_45, &InterceptionModel::_kpd_45);
        Internal(KPD_75, &InterceptionModel::_kpd_75);
        Internal(ALPHA, &InterceptionModel::_alpha);
        Internal(DOY, &InterceptionModel::_doy);

        //  external variables
        External(PAI, &InterceptionModel::_pai);
    }

    virtual ~InterceptionModel()
    {}

    void compute(double t, bool /* update */) {
        // parameters
        _Ta = _parameters.get(t).Temperature;
        _radiation = _parameters.get(t).Par;
        _doy = JulianCalculator::dayNumber(t);

        //Transform PAR in Global radiation
        _Rg = _radiation / _ec;

        //Transform day Rg in hourly Rg
        compute_HourlyRg();

        //Solar variables
        //Solar constant at the top of the atmosphere for a certain day
        _scd = 1370*(1+0.033*std::cos(360*_doy/365));
        //Seasonal offset of the solar height at a certain day
        _sinLD = std::sin(_latitudeRad)*std::sin(_declination);
        //Amplitude of sine of solar height at a certain day
        _cosLD = std::cos(_latitudeRad)*std::cos(_declination);
        //Day length
        _dayLength = 12+24/_pi*std::asin(_sinLD/_cosLD);
        //Integral solar height
        _solarHeight = 3600*(_dayLength*_sinLD+24/_pi*_cosLD*std::sqrt(1-std::pow(_sinLD/_cosLD,2)));
        //Daily extraterrestrial radiation
        _rExt = _scd * _solarHeight;
        //Atmospheric transmissivity
        _TransAtm = _Rg/(_rExt/1000000);

        //Compute diffuse and direct
        compute_diff_dir();

        //Compute canopy reflexion coefficient
        double rho_h = (1-std::sqrt(1-_zeta))/(1+std::sqrt(1-_zeta));


        //FOR EACH LAYER :
        _interc = 0;
        _Linterc = 0;
        for(int k = 0; k < _nbLayers; ++k) {
            //LAI of the current layer @TODO : compute LAI of layer
            _lai = _pai * (_density / 1.e4);

            //extinction coefficient for diffuse radiation
            _kpd_15 = compute_kd((15*_pi)/180);
            _kpd_45 = compute_kd((45*_pi)/180);
            _kpd_75 = compute_kd((75*_pi)/180);
            _alpha = std::sqrt(1-_zeta)*_lai;
            _kpd = -(1/_lai)*std::log(0.178*std::exp(-_kpd_15*_alpha)+0.514 * std::exp(-_kpd_45*_alpha)+0.308 * std::exp(-_kpd_75*_alpha));

            //compute diffuse PAr and direct PAR
            for(int i=0; i<24;++i) {
                _elevation_[i] = std::asin(_sinBeta_[i]);
                _kdr_bl_[i] = compute_kd(_elevation_[i]);
                _rho_cb_[i] = 1-std::exp((-2*rho_h*_kdr_bl_[i])/(1+_kdr_bl_[i]));
                //extinction coefficient for direct radiation
                _kpb_[i] = _kdr_bl_[i]*std::sqrt(1-_zeta);

                _diffusePar_[i] = _ec * (1- _rho_cd) * _rgHourly_diff_[i] * (1-std::exp(-(_kpd*_lai)));
                _directPar_[i] = _ec * (1- _rho_cb_[i]) * _rgHourly_dir_[i] * (1-std::exp(-(_kpb_[i]*_lai)));
                _interc_[i] = _diffusePar_[i] + _directPar_[i];
                _Linterc = _Linterc + _interc_[i];
            }
            _interc = _interc + _Linterc;
        }
    }


    //functions
    void compute_HourlyRg() {
        _declination = 23.45 * std::sin((2*_pi*(_doy+284))/366)*_pi/180;
        if(_latitudeRad+_declination < _pi/2 && _latitudeRad-_declination < _pi/2) { //normal day
            _ahs = std::acos(-std::tan(_latitudeRad)*std::tan(_declination));
            _sunrise = 12 - _ahs*12/_pi;
            _sunset = 12 + _ahs*12/_pi;
        } else if(_latitudeRad+_declination >= _pi/2) { //midnight sun
            _ahs = _pi;
            _sunrise = 0;
            _sunset = 24;
        } else if(_latitudeRad-_declination >= _pi/2) { //perpetual night
            _ahs = _pi;
            _sunrise = 12;
            _sunset = 12;
        }
        _dayDuration = _sunset - _sunrise;
        _mean = _sunrise + (_sunset-_sunrise)/2;
        _sigma = 0.15 * (_sunset-_sunrise);


        //@TODO : find truncated normal distribution code in c++
        const int nrolls=1000000;
        std::default_random_engine generator;
        generator.seed(1337);
        std::normal_distribution<double> distribution(_mean,_sigma);

        double p[24]={};
        for(int i=0; i<nrolls; ++i) {
            double number = distribution(generator);
            if(number>=0.0 && number < 24.0) {
                ++p[int(number)];
            }
        }
        for(int i=0;i<24;i++) {
            _rgHourly_[i] = (p[i]*_Rg)/nrolls;
            _sinBeta_[i] = std::max(0.0,std::sin(_latitudeRad)*std::sin(_declination)+std::cos(_latitudeRad)*std::cos(_declination)*cos(2*_pi*(i+12)/24));
            _rgHourly_diff_[i] = _propRgRd * _rgHourly_[i];
            _rgHourly_dir_[i] = _rgHourly_[i] - _rgHourly_diff_[i];
        }
    }

    void compute_diff_dir() {
        if (_TransAtm <= 0.07){
            _propRgRd = 1;
        } else if ((_TransAtm > 0.07) && (_TransAtm <= 0.35))  {
            _propRgRd = 1 - 2.3*std::pow(_TransAtm-0.07,2);
        } else if ((_TransAtm > 0.35) && (_TransAtm <= 0.75)) {
            _propRgRd = 1.33 - 1.46*_TransAtm;
        } else {
            _propRgRd = 0.23;
        }
        for(int i=0;i<24;i++) {
            _rgHourly_diff_[i] = _propRgRd * _rgHourly_[i];
            _rgHourly_dir_[i] = _rgHourly_[i] - _rgHourly_diff_[i];
        }
    }

    double compute_kd(double elevation) {
        double oav = 0;
        if(elevation>=(_leafAngle*_pi/180)) {
            oav = std::sin(elevation)*std::sin(_leafAngle*_pi/180);
        } else {
            oav = 2*((std::sin(elevation)*std::cos(_leafAngle*_pi/180) * std::asin(std::tan(elevation)/std::tan(_leafAngle*_pi/180)))+(std::sqrt(std::pow(std::sin(_leafAngle*_pi/180),2) - std::pow(std::sin(elevation),2))))/_pi;
        }
        if(sin(elevation) != 0) {
            return(oav/std::sin(elevation));
        } else {
            return(0);
        }
    }


    void init(double t, const ecomeristem::ModelParameters& parameters) {
        last_time = t-1;

        //parameters
        _parameters = parameters;
        _density = parameters.get("density");
        _ec = 0.48;
        _latitudeRad = 43.6167 * 3.141592653589793238462643383280/180;
        _pi = 3.141592653589793238462643383280;
        _zeta = 0.2;
        _leafAngle = 50;
        _rho_cd = 0.057;
        _nbLayers = 1;

        //  computed variables (internal)
        _Rg = 0;
        _declination = 0;
        _ahs = 0;
        _sunrise = 0;
        _sunset = 0;
        _dayDuration = 0;
        _mean = 0;
        _sigma = 0;
        _scd = 0;
        _sinLD = 0;
        _cosLD = 0;
        _dayLength = 0;
        _solarHeight = 0;
        _rExt = 0;
        _TransAtm = 0;
        _propRgRd = 1;
        _kpd = 0;
        _Linterc = 0;
        _interc = 0;
        _lai = 0;
        _kpd_15 = 0;
        _kpd_45 = 0;
        _kpd_75 = 0;
        _alpha = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;

    //  parameters
    double _ec;
    double _latitudeRad;
    double _pi;
    double _zeta;
    double _leafAngle;
    double _density;
    double _nbLayers;
    double _rho_cd;

    //  parameters(t)
    double _radiation;
    double _Ta;
    double _doy;

    //  internals - computed
    double _Rg;
    double _declination;
    double _ahs;
    double _sunrise;
    double _sunset;
    double _dayDuration;
    double _mean;
    double _sigma;

    double _rgHourly_[24];
    double _sinBeta_[24];
    double _elevation_[24];
    double _rgHourly_diff_[24];
    double _rgHourly_dir_[24];
    double _kdr_bl_[24];
    double _rho_cb_[24];
    double _kpb_[24];
    double _diffusePar_[24];
    double _directPar_[24];
    double _interc_[24];

    double _scd;
    double _sinLD;
    double _cosLD;
    double _dayLength;
    double _solarHeight;
    double _rExt;
    double _TransAtm;
    double _propRgRd;
    double _lai;
    double _kpd;
    double _interc;
    double _Linterc;
    double _kpd_15;
    double _kpd_45;
    double _kpd_75;
    double _alpha;


    //  externals
    double _pai;

};


} // namespace model
#endif
