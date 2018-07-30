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


#ifndef ASSIMILATION_MODEL_HPP
#define ASSIMILATION_MODEL_HPP

#include <defines.hpp>

namespace model {

class AssimilationModel : public AtomicModel < AssimilationModel >
{
public:
    enum internals { ASSIM, ASSIM_POT, INTERC, LAI, RESP_MAINT };

    enum externals { CSTR, FCSTR, FCSTRA, PAI, LEAFBIOMASS, INTERNODEBIOMASS, EXT_INTERC };


    AssimilationModel() {
        //  computed variables
        Internal(ASSIM, &AssimilationModel::_assim);
        Internal(ASSIM_POT, &AssimilationModel::_assim_pot);
        Internal(RESP_MAINT, &AssimilationModel::_resp_maint);
        Internal(INTERC, &AssimilationModel::_interc);
        Internal(LAI, &AssimilationModel::_lai);

        //  external variables
        External(CSTR, &AssimilationModel::_cstr);
        External(FCSTR, &AssimilationModel::_fcstr);
        External(FCSTRA, &AssimilationModel::_fcstrA);
        External(PAI, &AssimilationModel::_PAI);
        External(LEAFBIOMASS, &AssimilationModel::_LeafBiomass);
        External(INTERNODEBIOMASS, &AssimilationModel::_InternodeBiomass);
        External(EXT_INTERC, &AssimilationModel::_ext_interc);
    }

    virtual ~AssimilationModel()
    {}

    void compute(double t, bool /* update */) {
        // parameters
        _Ta = _parameters.get(t).Temperature;
        _radiation = _parameters.get(t).Par;

        //  lai
        _lai = _PAI * (_rolling_B + _rolling_A * _fcstr) * (_density / 1.e4);

        //  interc
        if(_intercmodel == 1) {
            _interc = 1. - std::exp(-_kdf * _lai);
        } else {
            _interc = _ext_interc;
        }

        //  assimPot
        if(_wbmodel == 1) {
            _assim_pot = std::pow(_cstr, _power_for_cstr) * _interc * _epsib * _radiation * _kpar;
        } else {
            _assim_pot = _fcstrA * _interc * _epsib * _radiation * _kpar;
        }

        //  respMaint
        _resp_maint = (_Kresp_leaf * _LeafBiomass + _Kresp_internode * _InternodeBiomass) * std::pow(2., (_Ta - _Tresp) / 10.);

        //  assim
        _assim = std::max(0., (_assim_pot / _density) - _resp_maint);
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        last_time = t-1;

        //parameters
        _parameters = parameters;

        //  parameters variables
        _density = parameters.get("density");
        _power_for_cstr = parameters.get("power_for_cstr");
        _kpar = 1 /* parameters.get("kpar") */;
        _epsib = parameters.get("Epsib");
        _kdf = parameters.get("Kdf");
        _rolling_A = parameters.get("Rolling_A");
        _rolling_B = parameters.get("Rolling_B");
        _Kresp_leaf = parameters.get("Kresp");
        _Kresp_internode = parameters.get("Kresp_internode");
        _Tresp = parameters.get("Tresp");
        _thresAssim = parameters.get("thresAssim");
        _wbmodel = parameters.get("wbmodel");
        _intercmodel = parameters.get("intercmodel");

        //  computed variables (internal)
        _assim = 0;
        _resp_maint = 0;
        _assim_pot = 0;
        _interc = 0;
        _lai = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;

    //  parameters
    double _density;
    double _power_for_cstr;
    double _kpar;
    double _epsib;
    double _kdf;
    double _rolling_A;
    double _rolling_B;
    double _Kresp_leaf;
    double _Kresp_internode;
    double _Tresp;
    double _intercmodel;

    double _thresAssim;
    double _wbmodel;

    //  parameters(t)
    double _radiation;
    double _Ta;

    //  internals - computed
    double _assim;
    double _resp_maint;
    double _assim_pot;
    double _interc;
    double _lai;

    //  externals
    double _cstr;
    double _fcstrA;
    double _fcstr;
    double _PAI;
    double _LeafBiomass;
    double _InternodeBiomass;
    double _ext_interc;
};

} // namespace model
#endif
