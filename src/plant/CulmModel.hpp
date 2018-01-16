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


#ifndef SAMPLE_ATOMIC_MODEL_HPP
#define SAMPLE_ATOMIC_MODEL_HPP

//#include <QDebug>
#include <defines.hpp>
#include <plant/processes/CulmStockModel.hpp>
#include <plant/processes/IctModel.hpp>
#include <plant/processes/ThermalTimeModelNG.hpp>
#include <plant/processes/ThermalTimeModel.hpp>
#include <plant/processes/CulmStockModelNG.hpp>
#include <plant/phytomer/PhytomerModel.hpp>
#include <plant/floralorgan/PanicleModel.hpp>
#include <plant/floralorgan/PeduncleModel.hpp>

namespace model {

class CulmModel : public CoupledModel < CulmModel >
{
public:

    enum submodels { PHYTOMERS, PANICLE, PEDUNCLE };


    enum internals { NB_LIG, NB_LIG_TOT, STEM_LEAF_PREDIM,
                     LEAF_BIOMASS_SUM, LEAF_LAST_DEMAND_SUM, LEAF_DEMAND_SUM,
                     INTERNODE_DEMAND_SUM, INTERNODE_LAST_DEMAND_SUM,
                     INTERNODE_BIOMASS_SUM, INTERNODE_LEN_SUM,
                     LEAF_BLADE_AREA_SUM,
                     REALLOC_BIOMASS_SUM, SENESC_DW_SUM,
                     LAST_LIGULATED_LEAF, LAST_LIGULATED_LEAF_LEN,
                     LAST_LEAF_BIOMASS_SUM, PANICLE_DAY_DEMAND, PEDUNCLE_BIOMASS,
                     PEDUNCLE_DAY_DEMAND, PEDUNCLE_LAST_DEMAND,
                     FIRST_LEAF_LEN, DELETED_SENESC_DW, PEDUNCLE_LEN, KILL_CULM,
                     LAST_LIGULATED_LEAF_BLADE_LEN, DELETED_REALLOC_BIOMASS, CULM_TEST_IC,
                     CULM_DEFICIT, CULM_STOCK, LAST_LEAF_INDEX, REALLOC_SUPPLY, PANICLE_WEIGHT,
                     LAST_LEAF_BLADE_AREA, NBLEAF,
                     DELETED_LEAF_NUMBER, LEAF_DELAY, SHEATH_LLL, REDUCTION_LER,
                     NB_INIT_LEAVES, NB_APP_LEAVES, NB_APP_LEAVES_TOT, BCPHY_CST, BCLIG_CST,
                     CULM_PHENO_STAGE, CULM_APP_STAGE, CULM_LIG_STAGE,
                     LIG_INDEX, MS_INDEX_CST, CULM_NBLEAF_PARAM2, CULM_PHENO_STAGE_CSTE,
                     PLASTO_INIT, PHYLLO_INIT, LIGULO_INIT, PLASTO, PHYLLO, LIGULO, TT_PLASTO,
                     TT_PHYLLO, TT_LIGULO };

    enum externals { PLANT_BOOL_CROSSED_PLASTO, DD, EDD, DELTA_T, FTSW, FCSTR, PLANT_PHENOSTAGE,
                     PLANT_APPSTAGE, PLANT_LIGSTAGE, PREDIM_LEAF_ON_MAINSTEM, SLA, PLANT_PHASE,
                     PLANT_STATE, TEST_IC, PLANT_STOCK, PLANT_DEFICIT, ASSIM, MGR,
                     LL_BL,IS_FIRST_DAY_PI, MS_PHYT_INDEX, MS_SHEATH_LLL };


    CulmModel(int index):
        _index(index), _is_first_culm(index == 1),
        _culm_stock_model(new CulmStockModelNG),
        _culm_thermaltime_model(new ThermalTimeModel),
        _culm_ictmodel(new IctModel),
        _culm_thermaltime_modelNG(new ThermalTimeModelNG)
    {

        Internal(NB_LIG, &CulmModel::_nb_lig);
        Internal(STEM_LEAF_PREDIM, &CulmModel::_stem_leaf_predim);
        Internal(LEAF_BIOMASS_SUM, &CulmModel::_leaf_biomass_sum);
        Internal(LEAF_LAST_DEMAND_SUM, &CulmModel::_leaf_last_demand_sum);
        Internal(LEAF_DEMAND_SUM, &CulmModel::_leaf_demand_sum);
        Internal(INTERNODE_LAST_DEMAND_SUM, &CulmModel::_internode_last_demand_sum);
        Internal(INTERNODE_DEMAND_SUM, &CulmModel::_internode_demand_sum);
        Internal(INTERNODE_BIOMASS_SUM, &CulmModel::_internode_biomass_sum);
        Internal(INTERNODE_LEN_SUM, &CulmModel::_internode_len_sum);
        Internal(LEAF_BLADE_AREA_SUM, &CulmModel::_leaf_blade_area_sum);
        Internal(REALLOC_BIOMASS_SUM, &CulmModel::_realloc_biomass_sum);
        Internal(SENESC_DW_SUM, &CulmModel::_senesc_dw_sum);
        Internal(LAST_LIGULATED_LEAF, &CulmModel::_last_ligulated_leaf);
        Internal(LAST_LIGULATED_LEAF_LEN, &CulmModel::_last_ligulated_leaf_len);
        Internal(LAST_LEAF_BIOMASS_SUM, &CulmModel::_last_leaf_biomass_sum);
        Internal(PANICLE_DAY_DEMAND, &CulmModel::_panicle_day_demand);
        Internal(PEDUNCLE_BIOMASS, &CulmModel::_peduncle_biomass);
        Internal(PEDUNCLE_DAY_DEMAND, &CulmModel::_peduncle_day_demand);
        Internal(PEDUNCLE_LAST_DEMAND, &CulmModel::_peduncle_last_demand);
        Internal(FIRST_LEAF_LEN, &CulmModel::_first_leaf_len);
        Internal(DELETED_SENESC_DW, &CulmModel::_deleted_senesc_dw);
        Internal(NB_LIG_TOT, &CulmModel::_nb_lig_tot);
        Internal(PEDUNCLE_LEN, &CulmModel::_peduncle_len);
        Internal(KILL_CULM, &CulmModel::_kill_culm);
        Internal(LAST_LIGULATED_LEAF_BLADE_LEN, &CulmModel::_last_ligulated_leaf_blade_len);
        Internal(DELETED_REALLOC_BIOMASS, &CulmModel::_deleted_realloc_biomass);
        Internal(CULM_TEST_IC, &CulmModel::_culm_test_ic);
        Internal(CULM_DEFICIT, &CulmModel::_culm_deficit);
        Internal(CULM_STOCK, &CulmModel::_culm_stock);
        Internal(LAST_LEAF_INDEX, &CulmModel::_last_leaf_index);
        Internal(REALLOC_SUPPLY, &CulmModel::_realloc_supply);
        Internal(PANICLE_WEIGHT, &CulmModel::_panicle_weight);
        Internal(LAST_LEAF_BLADE_AREA, &CulmModel::_last_leaf_blade_area);
        Internal(NBLEAF, &CulmModel::_nb_leaf);
        Internal(DELETED_LEAF_NUMBER, &CulmModel::_deleted_leaf_number);
        Internal(LEAF_DELAY, &CulmModel::_leaf_delay);
        Internal(SHEATH_LLL, &CulmModel::_sheath_LLL);
        Internal(REDUCTION_LER, &CulmModel::_reductionLER);
        Internal(NB_INIT_LEAVES, &CulmModel::_nb_init_leaves);
        Internal(NB_APP_LEAVES, &CulmModel::_nb_app_leaves);
        Internal(NB_APP_LEAVES_TOT, &CulmModel::_nb_app_leaves_tot);
        Internal(CULM_PHENO_STAGE, &CulmModel::_culm_phenostage);
        Internal(CULM_APP_STAGE, &CulmModel::_culm_appstage);
        Internal(CULM_LIG_STAGE, &CulmModel::_culm_ligstage);
        Internal(LIG_INDEX, &CulmModel::_lig_index);
        Internal(CULM_NBLEAF_PARAM2, &CulmModel::_culm_nbleaf_param_2);
        Internal(CULM_PHENO_STAGE_CSTE, &CulmModel::_culm_phenostage_cste);
        Internal(PLASTO_INIT, &CulmModel::_plasto_init);
        Internal(PHYLLO_INIT, &CulmModel::_phyllo_init);
        Internal(LIGULO_INIT, &CulmModel::_ligulo_init);
        Internal(PLASTO, &CulmModel::_plasto);
        Internal(PHYLLO, &CulmModel::_phyllo);
        Internal(LIGULO, &CulmModel::_ligulo);
        Internal(TT_PLASTO, &CulmModel::_tt_plasto);
        Internal(TT_PHYLLO, &CulmModel::_tt_phyllo);
        Internal(TT_LIGULO, &CulmModel::_tt_ligulo);

        //    externals
        External(PLANT_BOOL_CROSSED_PLASTO, &CulmModel::_bool_crossed_plasto);
        External(DD, &CulmModel::_dd);
        External(EDD, &CulmModel::_edd);
        External(DELTA_T, &CulmModel::_delta_t);
        External(FTSW, &CulmModel::_ftsw);
        External(FCSTR, &CulmModel::_fcstr);
        External(PLANT_PHENOSTAGE, &CulmModel::_plant_phenostage);
        External(PLANT_APPSTAGE, &CulmModel::_plant_appstage);
        External(PLANT_LIGSTAGE, &CulmModel::_plant_ligstage);
        External(PREDIM_LEAF_ON_MAINSTEM, &CulmModel::_predim_leaf_on_mainstem);
        External(SLA, &CulmModel::_sla);
        External(PLANT_STATE, &CulmModel::_plant_state);
        External(PLANT_PHASE, &CulmModel::_plant_phase);
        External(TEST_IC, &CulmModel::_test_ic);
        External(PLANT_STOCK, &CulmModel::_plant_stock);
        External(PLANT_DEFICIT, &CulmModel::_plant_deficit);
        External(ASSIM, &CulmModel::_assim);
        External(MGR, &CulmModel::_MGR);
        External(LL_BL, &CulmModel::_LL_BL);
        External(IS_FIRST_DAY_PI, &CulmModel::_is_first_day_pi);
        External(MS_SHEATH_LLL, &CulmModel::_ms_sheath_LLL);
        External(MS_PHYT_INDEX, &CulmModel::_ms_phyt_index);
    }

    virtual ~CulmModel()
    {
        _culm_stock_model.reset(nullptr);
        _culm_ictmodel.reset(nullptr);
        _culm_thermaltime_modelNG.reset(nullptr);
        _culm_thermaltime_model.reset(nullptr);
        _panicle_model.reset(nullptr);
        _peduncle_model.reset(nullptr);

        auto it = _phytomer_models.begin();
        while (it != _phytomer_models.end()) {
            delete *it;
            ++it;
        }
    }


    bool is_phytomer_creatable() {
        return (_culm_phase == culm::VEGETATIVE
                || _culm_phase == culm::ELONG
                || _culm_phase == culm::PI
                )
                && (get_phytomer_number() < _maxleaves)
                && (_culm_deficit + _culm_stock >= 0);
    }

    void step_state(double t) {

        if(_plant_phase == plant::PI && _lag == false && _is_lagged == false) {
            _lag = true;
            _is_lagged = true;
            if(_culm_phenostage_at_lag == 0) {
                _culm_phenostage_at_lag = _culm_phenostage;
            }
        }

        if(_lag) {
            if(_culm_phenostage == _culm_phenostage_at_lag + _coeff_pi_lag && get_phytomer_number() >= 3) {
                _culm_phase = culm::PI;
                if(_culm_phenostage_at_pi == 0) {
                    _culm_phenostage_at_pi = _culm_phenostage;
                }
                _lag = false;
            }
        }

        if( _plant_phase == plant::PI && _is_first_culm) {
            _culm_phase = culm::PI;
            if(_culm_phenostage_at_pi == 0) {
                _culm_phenostage_at_pi = _culm_phenostage;
            }
        }

        switch( _culm_phase  ) {
        case culm::INITIAL: {
            _culm_phase  = culm::VEGETATIVE;
            break;
        }
        case culm::VEGETATIVE: {
            if(_plant_phase == plant::ELONG) {
                if(_nb_lig > 0) {
                    _culm_phase  = culm::ELONG;
                }
            }
            break;
        }
        case culm::ELONG: {
            break;
        }
        case culm::PI: {
            _last_phase = _culm_phase ;
            if( _plant_state & plant::NEW_PHYTOMER_AVAILABLE) {
                if( _plant_phase == plant::PI) {
                    if(!_started_PI) {
                        _panicle_model = std::unique_ptr<PanicleModel>(new PanicleModel());
                        setsubmodel(PANICLE, _panicle_model.get());
                        _panicle_model->init(t, _parameters);
                        _started_PI = true;
                    }
                }
            }
            if (_culm_ligstage == _maxleaves + 1) {
                if(_is_first_culm and _plant_phase == plant::PI) {
                    return;
                }
                _peduncle_model = std::unique_ptr<PeduncleModel>(new PeduncleModel(_index, _is_first_culm));
                setsubmodel(PEDUNCLE, _peduncle_model.get());
                _peduncle_model->init(t, _parameters);
                _culm_phase = culm::PRE_FLO;
                _culm_phenostage_at_pre_flo = _culm_phenostage;
            }
            break;
        }
        case culm::PRE_FLO: {
            _last_phase = _culm_phase ;
            if(_culm_ligstage == _maxleaves + 4) { //+1 to pre_flo; +3 to end peduncle growth)
                _culm_phase = culm::FLO;
            }
            break;
        }
        }
    }

    void compute(double t, bool /* update */) {
        std::string date = artis::utils::DateTime::toJulianDayFmt(t, artis::utils::DATE_FORMAT_YMD);
        if(_kill_culm or _plant_phase == plant::DEAD) {
            _nb_lig = 0;
            _nb_lig_tot = 0;
            _leaf_biomass_sum = 0;
            _last_leaf_biomass_sum = 0;
            _leaf_last_demand_sum = 0;
            _leaf_demand_sum = 0;
            _internode_last_demand_sum = 0;
            _internode_demand_sum = 0;
            _internode_biomass_sum = 0;
            _internode_len_sum = 0;
            _leaf_blade_area_sum = 0;
            _last_ligulated_leaf = -1;
            _last_ligulated_leaf_len = 0;
            _realloc_biomass_sum = 0;
            _deleted_realloc_biomass = 0;
            _realloc_supply = 0;
            _nb_leaf = 0;
            _culm_phase = culm::DEAD;
            _leaf_delay = 0;
            _kill_culm = true;
            auto it = _phytomer_models.begin();
            while (it != _phytomer_models.end()) {
                if(!(*it)->is_leaf_dead(t-1)) {
                    (*it)->kill_leaf(t);
                }
                ++it;
            }
            return;
        }

        //thermaltimevalues update
        if(_plant_phenostage >= _nb_leaf_param2 - 1 and _plant_stock >= 0) {
            _plasto = _plasto_init * _coeff_Plasto_PI;
            _phyllo = _phyllo_init * _coeff_Phyllo_PI;
            _ligulo = _ligulo_init * _coeff_Ligulo_PI;
            _tt_plasto = _plasto;
        }
        if(_plant_appstage >= _nb_leaf_param2 - 1) {
            _tt_phyllo = _phyllo;
        }
        if(_plant_ligstage >= _nb_leaf_param2 - 1) {
            _tt_ligulo = _ligulo;
        }

        if(_creation_date == t) {
            _culm_thermaltime_model->put(t, ThermalTimeModel::DELTA_T, _delta_t);
            _culm_thermaltime_model->put(t, ThermalTimeModel::PLASTO, _plasto);
            _culm_thermaltime_model->put(t, ThermalTimeModel::PHYLLO, _phyllo);
            _culm_thermaltime_model->put(t, ThermalTimeModel::LIGULO, _ligulo);
            _culm_thermaltime_model->put(t, ThermalTimeModel::PLASTO_DELAY, _leaf_delay);
            _culm_thermaltime_model->put(t, ThermalTimeModel::STOCK, _plant_stock);
            compute_thermaltime(t);
        }
        if(_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE or _is_first_day_pi) {
            if(_culm_thermaltime_model) {
                _culm_DD = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::DD);
                _culm_EDD = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::EDD);
                _culm_bool_crossed_plasto = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::BOOL_CROSSED_PLASTO);
                _culm_bool_crossed_phyllo = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::BOOL_CROSSED_PHYLLO);
                _culm_bool_crossed_ligulo = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::BOOL_CROSSED_LIGULO);
                _culm_plasto_visu = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::PHYLLO_VISU);
                _culm_ligulo_visu = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::LIGULO_VISU);
                _culm_phyllo_visu = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::PHYLLO_VISU);
                _culm_DD_phyllo = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::DD_PHYLLO);
                _culm_DD_ligulo = _culm_thermaltime_model->get < double >(t, ThermalTimeModel::DD_LIGULO);
                _culm_phenostage = _culm_thermaltime_model->get < int >(t, ThermalTimeModel::PHENO_STAGE);
                _culm_appstage = _culm_thermaltime_model->get < int >(t, ThermalTimeModel::APP_STAGE);
                _culm_ligstage = _culm_thermaltime_model->get < int >(t, ThermalTimeModel::LIG_STAGE);
                _culm_phenostage_cste = _culm_phenostage;
                _culm_appstage_cste = _culm_appstage;
                _culm_ligstage_cste = _culm_ligstage;
            }
        } else {
            _culm_DD = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_DD);
            _culm_EDD = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::EDD);
            _culm_bool_crossed_plasto = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PLASTO);
            _culm_bool_crossed_phyllo = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PHYLLO);
            _culm_bool_crossed_ligulo = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_LIGULO);
            _culm_plasto_visu = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_PHYLLO_VISU);
            _culm_ligulo_visu = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_LIGULO_VISU);
            _culm_phyllo_visu = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_PHYLLO_VISU);
            _culm_DD_phyllo = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_DD_PHYLLO);
            _culm_DD_ligulo = _culm_thermaltime_modelNG->get < double >(t, ThermalTimeModelNG::CULM_DD_LIGULO);
            _culm_phenostage = _culm_phenostage_cste + _culm_thermaltime_modelNG->get < int >(t, ThermalTimeModelNG::PHENO_STAGE);
            _culm_appstage = _culm_appstage_cste + _culm_thermaltime_modelNG->get < int >(t, ThermalTimeModelNG::APP_STAGE);
            _culm_ligstage = _culm_ligstage_cste + _culm_thermaltime_modelNG->get < int >(t, ThermalTimeModelNG::LIG_STAGE);
        }

        //phytomer creation
        if(_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE) {
            if( ( _plant_state & plant::NEW_PHYTOMER_AVAILABLE ) && is_phytomer_creatable()) {
                create_phytomer(t);
            }
        } else {
            if(is_phytomer_creatable() and _culm_bool_crossed_plasto >= 0) {
                create_phytomer(t);
            }
        }

        //last index before elongation
        if(_culm_phase == culm::VEGETATIVE or _culm_phase == culm::INITIAL) {
            _last_leaf_index = _culm_ligstage;
        }

        step_state(t);

        //nb leaf param2 for specific culm
        if(_plant_phenostage == _nb_leaf_param2 and _bool_crossed_plasto >= 0 and _plant_stock > 0) {
            _culm_nbleaf_param_2 = _phytomer_models.size();
        }

        //kill leaf
        _realloc_biomass_sum = 0;
        _deleted_senesc_dw = 0;
        _deleted_realloc_biomass = 0;
        _realloc_supply = 0;
        if(_plant_phase != plant::INITIAL and _plant_phase != plant::VEGETATIVE) {
            if(_is_first_day_pi or t == _creation_date) {
                if((_plant_deficit * (_leaf_biomass_sum / _plant_leaf_biomass_sum)) + (_plant_stock * (_leaf_biomass_sum / _plant_leaf_biomass_sum)) < 0) {
                    delete_leafNG(t, get_first_alive_leaf_index2(t));
                    _realloc_biomass_sum += _deleted_realloc_biomass;
                }
            } else if(! _kill_culm) {
                if((_culm_stock_model->get < double >(t-1, CulmStockModelNG::CULM_STOCK)) + (_culm_stock_model->get < double >(t-1, CulmStockModelNG::CULM_DEFICIT)) < 0) {
                    delete_leafNG(t, get_first_alive_leaf_index2(t));
                    _realloc_biomass_sum += _deleted_realloc_biomass;
                }
            }
        }

        auto it = _phytomer_models.begin();
        std::deque < PhytomerModel* >::iterator previous_it;
        int i = 0;
        _nb_lig = 0;
        _nb_lig_tot = 0;
        _leaf_biomass_sum = 0;
        _last_leaf_biomass_sum = 0;
        _leaf_last_demand_sum = 0;
        _leaf_demand_sum = 0;
        _internode_last_demand_sum = 0;
        _internode_demand_sum = 0;
        _internode_biomass_sum = 0;
        _internode_len_sum = 0;
        _leaf_blade_area_sum = 0;
        _last_ligulated_leaf = -1;
        _last_ligulated_leaf_len = 0;
        _nb_leaf = 0;
        _nb_init_leaves = 0;
        _nb_app_leaves = 0;
        _nb_app_leaves_tot = 0;
        _reductionLER = 1.;
        _sheath_LLL = _ms_sheath_LLL;

        while (it != _phytomer_models.end()) {
            //Phytomers
            compute_phytomers(it, previous_it, i, t);
            //Sum
            compute_vars(it, previous_it, i, t);
            //GetLastINnonVegetative
            get_nonvegetative_in(it, t);

            if(i == 0 and not (*it)->is_leaf_dead(t)) {
                _first_leaf_len = (*it)->leaf()->get < double >(t, LeafModel::BLADE_LEN);
                if ((*it)->is_leaf_lig(t) and t == (*it)->leaf()->get < double >(t, LeafModel::LIG_T)) {
                    _last_ligulated_leaf = i;
                    _last_ligulated_leaf_len = (*it)->leaf()->get < double >(t, LeafModel::LEAF_LEN);
                    _last_ligulated_leaf_blade_len = (*it)->leaf()->get < double >(t, LeafModel::BLADE_LEN);
                }
            }

            previous_it = it;
            ++it;
            ++i;
        }

        //Floral_organs
        if(_panicle_model.get()) {
            _panicle_model->put (t, PanicleModel::DELTA_T, _delta_t);
            _panicle_model->put < plant::plant_phase >(t, PanicleModel::PLANT_PHASE, _plant_phase);
            _panicle_model->put(t, PanicleModel::FCSTR, _fcstr);
            _panicle_model->put(t, PanicleModel::TEST_IC, _test_ic);
            (*_panicle_model)(t);
            _panicle_day_demand = _panicle_model->get < double >(t, PanicleModel::DAY_DEMAND);
            _panicle_weight = _panicle_model->get < double >(t, PanicleModel::WEIGHT);
        }

        if(_peduncle_model.get()) {
            _peduncle_model->put < plant::plant_phase >(t, PeduncleModel::PLANT_PHASE, _plant_phase);
            _peduncle_model->put < culm::culm_phase >(t, PeduncleModel::CULM_PHASE, _culm_phase);
            _peduncle_model->put (t, PeduncleModel::INTER_PREDIM, _peduncle_inerlen_predim);
            _peduncle_model->put (t, PeduncleModel::INTER_DIAM, _peduncle_inerdiam_predim);
            _peduncle_model->put(t, PeduncleModel::FTSW, _ftsw);
            _peduncle_model->put(t, PeduncleModel::EDD, _culm_EDD);
            _peduncle_model->put (t, PeduncleModel::DELTA_T, _delta_t);
            _peduncle_model->put (t, PeduncleModel::PLASTO, _plasto);
            _peduncle_model->put (t, PeduncleModel::LIGULO, _ligulo);
            _peduncle_model->put (t, PeduncleModel::FCSTR, _fcstr);
            _peduncle_model->put (t, PeduncleModel::TEST_IC, _test_ic);
            (*_peduncle_model)(t);
        }

        if(_peduncle_model.get()) {
            _peduncle_last_demand = _peduncle_model->get < double >(t, PeduncleModel::LAST_DEMAND);
            _peduncle_day_demand = _peduncle_model->get < double >(t, PeduncleModel::DEMAND);
            _peduncle_len = _peduncle_model->get < double >(t, PeduncleModel::LENGTH);
            _peduncle_biomass = _peduncle_model->get < double >(t, PeduncleModel::BIOMASS);
        }

        //Plastodelay
        _leaf_delay = _delta_t * (-1. + _reductionLER);

        //kill culms without leaves
        if(get_alive_phytomer_number() <= 0) {
            _kill_culm = true;
        }
    }

    void get_nonvegetative_in(std::deque < PhytomerModel* >::iterator it, double t) {
        if(_peduncle_model.get()) {
            if(((*it)->internode()->get < double >(t, InternodeModel::INTERNODE_LEN)) > 0) {
                _peduncle_inerlen_predim = (*it)->internode()->get<double>(t, InternodeModel::INTERNODE_PREDIM);
                _peduncle_inerdiam_predim = (*it)->internode()->get<double>(t, InternodeModel::INTER_DIAMETER);
            }
        }
    }

    void compute_stock(double t) {
        _culm_stock_model->put(t, CulmStockModelNG::PLANT_STOCK, _plant_stock);
        _culm_stock_model->put(t, CulmStockModelNG::LEAF_BIOMASS_SUM, _leaf_biomass_sum);
        _culm_stock_model->put(t, CulmStockModelNG::INTERNODE_BIOMASS_SUM, _internode_biomass_sum);
        _culm_stock_model->put(t, CulmStockModelNG::LEAF_DEMAND_SUM, _leaf_demand_sum);
        _culm_stock_model->put(t, CulmStockModelNG::INTERNODE_DEMAND_SUM, _internode_demand_sum);
        _culm_stock_model->put(t, CulmStockModelNG::LAST_DEMAND, _leaf_last_demand_sum + _internode_last_demand_sum + _peduncle_last_demand);
        _culm_stock_model->put(t, CulmStockModelNG::REALLOC_BIOMASS, _realloc_biomass_sum);
        _culm_stock_model->put(t, CulmStockModelNG::PLANT_PHASE, _plant_phase);
        _culm_stock_model->put(t, CulmStockModelNG::PANICLE_DEMAND, _panicle_day_demand);
        _culm_stock_model->put(t, CulmStockModelNG::IS_FIRST_DAY_OF_INDIVIDUALIZATION, _is_first_day_pi);
        _culm_stock_model->put(t, CulmStockModelNG::PEDUNCLE_DEMAND, _peduncle_day_demand);
        _culm_stock_model->put(t, CulmStockModelNG::KILL_CULM, _kill_culm);
        (*_culm_stock_model)(t);
        _culm_test_ic = std::min(1.,sqrt(_culm_stock_model->get < double >(t, CulmStockModelNG::CULM_IC)));
        _culm_deficit = _culm_stock_model->get < double >(t, CulmStockModelNG::CULM_DEFICIT);
        _culm_stock = _culm_stock_model->get < double >(t, CulmStockModelNG::CULM_STOCK);
    }

    void compute_ictmodel(double t) {
        _culm_ictmodel->put(t, IctModel::BOOL_CROSSED_PHYLLO, _culm_bool_crossed_phyllo);
        (*_culm_ictmodel)(t);
    }

    void compute_thermaltime(double t) {
        _culm_thermaltime_model->put(t, ThermalTimeModel::IS_FIRST_CULM, _is_first_culm);
        _culm_thermaltime_model->put(t, ThermalTimeModel::PLASTO, _tt_plasto);
        _culm_thermaltime_model->put(t, ThermalTimeModel::PHYLLO, _tt_phyllo);
        _culm_thermaltime_model->put(t, ThermalTimeModel::LIGULO, _tt_ligulo);
        (*_culm_thermaltime_model)(t);
    }

    void compute_thermaltimeNG(double t) {
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::CULM_STOCK, _culm_stock);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::CULM_DEFICIT, _culm_deficit);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::PLASTO, _tt_plasto);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::PHYLLO, _tt_phyllo);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::LIGULO, _tt_ligulo);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::IS_FIRST_DAY_OF_INDIVIDUALIZATION, _is_first_day_pi);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::PLASTO_DELAY, _leaf_delay);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::PLASTO_VISU, _culm_plasto_visu);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::LIGULO_VISU, _culm_ligulo_visu);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::PHYLLO_VISU, _culm_phyllo_visu);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::DD, _culm_DD);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::DD_PHYLLO, _culm_DD_phyllo);
        _culm_thermaltime_modelNG->put(t, ThermalTimeModelNG::DD_LIGULO, _culm_DD_ligulo);
        (*_culm_thermaltime_modelNG)(t);
    }

    void compute_phytomers(std::deque < PhytomerModel* >::iterator it, std::deque < PhytomerModel* >::iterator previous_it, int i, double t) {
        (*it)->put(t, PhytomerModel::DD, _culm_DD);
        (*it)->put(t, PhytomerModel::DELTA_T, _delta_t);
        (*it)->put(t, PhytomerModel::FTSW, _ftsw);
        (*it)->put(t, PhytomerModel::FCSTR, _fcstr);
        (*it)->put(t, PhytomerModel::SLA, _sla);
        (*it)->put < plant::plant_state >(t, PhytomerModel::PLANT_STATE, _plant_state);
        (*it)->put < plant::plant_phase >(t, PhytomerModel::PLANT_PHASE, _plant_phase);
        if(_plant_phase != plant::INITIAL and _plant_phase != plant::VEGETATIVE and !(_is_first_day_pi)) {
            (*it)->put(t, PhytomerModel::TEST_IC, _culm_test_ic);
        } else {
            (*it)->put(t, PhytomerModel::TEST_IC, _test_ic);
        }
        (*it)->leaf()->put(t, LeafModel::MGR, _MGR);
        (*it)->leaf()->put(t, LeafModel::CULM_DEFICIT, _culm_deficit);
        (*it)->leaf()->put(t, LeafModel::CULM_STOCK, _culm_stock);
        (*it)->leaf()->put(t, LeafModel::SHEATH_LLL, _sheath_LLL);
        (*it)->leaf()->put(t, LeafModel::CULM_NBLEAF_PARAM2, _culm_nbleaf_param_2);
        (*it)->internode()->put(t, InternodeModel::CULM_DEFICIT, _culm_deficit);
        (*it)->internode()->put(t, InternodeModel::CULM_STOCK, _culm_stock);
        (*it)->internode()->put(t, InternodeModel::CULM_PHASE, _culm_phase);
        (*it)->internode()->put(t, InternodeModel::PLASTO, _plasto);
        (*it)->internode()->put(t, InternodeModel::LIGULO, _ligulo);
        (*it)->internode()->put(t, InternodeModel::NB_LIG, _nb_lig);
        (*it)->internode()->put(t, InternodeModel::BOOL_CROSSED_PLASTO, _culm_bool_crossed_plasto);
        (*it)->internode()->put(t, InternodeModel::LAST_LEAF_INDEX, _last_leaf_index);
        (*it)->internode()->put(t, InternodeModel::PHENOSTAGE, _plant_phenostage);
        (*it)->internode()->put(t, InternodeModel::LIGSTAGE, _plant_ligstage);
        if(i == 0) {
            (*it)->internode()->put(t, InternodeModel::PREVIOUS_IN_PREDIM, 0.);
        } else {
            (*it)->internode()->put(t, InternodeModel::PREVIOUS_IN_PREDIM, (*previous_it)->internode()->get < double >(t, InternodeModel::INTERNODE_PREDIM));
        }
        if (_is_first_culm) {
            if (i == 0) {
                (*it)->put(t, PhytomerModel::PREDIM_LEAF_ON_MAINSTEM, 0.);
            } else {
                (*it)->put(t, PhytomerModel::PREDIM_LEAF_ON_MAINSTEM,
                           (*previous_it)->get < double, LeafModel >(t, PhytomerModel::LEAF_PREDIM));
            }
        } else {
            (*it)->put(t, PhytomerModel::PREDIM_LEAF_ON_MAINSTEM, _predim_leaf_on_mainstem);
        }

        if (i == 0) {
            (*it)->put(t, PhytomerModel::PREDIM_PREVIOUS_LEAF, 0.);
        } else {
            (*it)->put(t, PhytomerModel::PREDIM_PREVIOUS_LEAF,
                       (*previous_it)->get < double, LeafModel >(t, PhytomerModel::LEAF_PREDIM));
        }
        (**it)(t);
    }

    void compute_vars(std::deque < PhytomerModel* >::iterator it, std::deque < PhytomerModel* >::iterator previous_it, int i, double t) {
        if((*it)->is_leaf_lig(t)) {
            _lig_index = i+1;
        }
        if (!(*it)->is_leaf_dead(t)) {
            if((*it)->is_leaf_lig(t)) {
                ++ _nb_lig;
                ++ _nb_lig_tot;
                double ratio = (*it)->leaf()->get < double >(t, LeafModel::BLADE_AREA) / (*it)->leaf()->get < double >(t, LeafModel::LAST_BLADE_AREA);
                if(ratio >= 0.5) {
                    _nb_leaf = _nb_leaf + 1;
                }
                ++ _nb_app_leaves_tot;
            } else {
                _nb_leaf = _nb_leaf + 1;
                _reductionLER = (*it)->leaf()->get < double >(t, LeafModel::REDUCTION_LER);
                if((*it)->is_leaf_app(t)) {
                    _nb_app_leaves = _nb_app_leaves + 1;
                    ++ _nb_app_leaves_tot;
                } else {
                    _nb_init_leaves = _nb_init_leaves + 1;
                }
            }

            if (_index == 1) {
                _stem_leaf_predim = (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_PREDIM);
                if(!((*it)->is_leaf_dead(t)) and (*it)->is_leaf_lig(t)) {
                    _last_leaf_blade_area = (*it)->leaf()->get < double >(t, LeafModel::LAST_BLADE_AREA);
                }
            }
            if (_index == 1 and (*it)->is_leaf_lig(t) and t != (*it)->leaf()->get < double >(t, LeafModel::LIG_T)) {
                _last_ligulated_leaf = i;
                _last_ligulated_leaf_len = (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_LEN);
                _last_ligulated_leaf_blade_len = (*it)->leaf()->get < double >(t, LeafModel::BLADE_LEN);
            }
        } else {
            ++_nb_lig_tot;
            if((*it)->is_leaf_ligged(t) or (*it)->is_leaf_apped(t)) {
                ++ _nb_app_leaves_tot;
            }
        }
        if(i == 0) {
            _sheath_LLL = ((*it)->leaf()->get < double >(t, LeafModel::LEAF_PREDIM)) - ((1 - (1 / _LL_BL))*(*it)->leaf()->get < double >(t, LeafModel::LEAF_PREDIM));
        }
        if((*it)->is_leaf_lig(t)) {
            _sheath_LLL = (*it)->leaf()->get < double >(t, LeafModel::LEAF_LEN) - (*it)->leaf()->get < double >(t, LeafModel::BLADE_LEN);
        }
        _leaf_biomass_sum += (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_BIOMASS);

        if ((*it)->leaf()->get < double >(t, LeafModel::LAST_LEAF_BIOMASS) == 0) {
            _last_leaf_biomass_sum += (*it)->leaf()->get < double >(t, LeafModel::BIOMASS);
        } else {
            _last_leaf_biomass_sum += (*it)->leaf()->get < double >(t, LeafModel::LAST_LEAF_BIOMASS);
        }
        _leaf_last_demand_sum +=
                (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_LAST_DEMAND);
        _leaf_demand_sum += (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_DEMAND);
        _internode_last_demand_sum +=
                (*it)->get < double, InternodeModel >(t, PhytomerModel::INTERNODE_LAST_DEMAND);
        _internode_demand_sum +=
                (*it)->get < double, InternodeModel >(
                    t, PhytomerModel::INTERNODE_DEMAND);
        _internode_biomass_sum += (*it)->get < double, InternodeModel >(
                    t, PhytomerModel::INTERNODE_BIOMASS);
        _internode_len_sum += (*it)->get < double, InternodeModel >(
                    t, PhytomerModel::INTERNODE_LEN);
        _leaf_blade_area_sum +=
                (*it)->get < double, LeafModel >(
                    t, PhytomerModel::LEAF_BLADE_AREA);

        _realloc_biomass_sum +=
                (*it)->get < double, LeafModel >(
                    t, PhytomerModel::REALLOC_BIOMASS);
        _senesc_dw_sum +=
                (*it)->get < double, LeafModel >(
                    t, PhytomerModel::SENESC_DW);
    }

    void create_phytomer(double t)
    {
        if (t != _parameters.beginDate) {
            int index;
            if (_phytomer_models.empty()) {
                index = 1;
            } else {
                index = _phytomer_models.back()->get_index() + 1;
            }

            PhytomerModel* phytomer = new PhytomerModel(index, _is_first_culm, _plasto, _phyllo, _ligulo, _LL_BL);
            setsubmodel(PHYTOMERS, phytomer);
            phytomer->init(t, _parameters);
            _phytomer_models.push_back(phytomer);
        }
    }

    int get_phytomer_number() const
    { return _phytomer_models.size(); }

    int get_app_phytomer_number(double t) const {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        int i = 0;
        while (it != _phytomer_models.end()) {
            if ((*it)->is_leaf_dead(t-1) or (*it)->is_leaf_lig(t-1) or (*it)->is_leaf_app(t-1)) {
                ++i;
            }
            ++it;
        }
        return i;
    }

    int get_alive_phytomer_number() const
    { return _phytomer_models.size() - _deleted_leaf_number; }

    int get_dead_phytomer_number(double t) const {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        int i = 0;
        while (it != _phytomer_models.end()) {
            if ((*it)->is_leaf_dead(t)) {
                ++i;
            }
            ++it;
        }
        return i;
    }

    CulmStockModelNG * stock_model() const {
        return _culm_stock_model.get();
    }

    IctModel * ictmodel() const {
            return _culm_ictmodel.get();
    }

    ThermalTimeModel * thermaltime_model() const {
        return _culm_thermaltime_model.get();
    }

    ThermalTimeModelNG * thermaltime_modelNG() const {
        return _culm_thermaltime_modelNG.get();
    }

    void delete_leaf(double t, int index, double leaf_biomass, double internode_biomass)
    {
        _deleted_senesc_dw += (1 - _realocationCoeff) * leaf_biomass;
        _phytomer_models[index]->kill_leaf(t);
        ++_deleted_leaf_number;
        if (get_alive_phytomer_number() == 0) {
            _kill_culm = true;
            _deleted_senesc_dw += (1 - _realocationCoeff) * internode_biomass;
        }
    }

    void delete_leafNG(double t, int index)
    {
        if(index != -1) {
            _deleted_senesc_dw = (1 - _realocationCoeff) * _phytomer_models[index]->leaf()->get < double >(t-1, LeafModel::BIOMASS);
            _deleted_realloc_biomass = (_realocationCoeff) * _phytomer_models[index]->leaf()->get < double >(t-1, LeafModel::BIOMASS);
            _phytomer_models[index]->kill_leaf(t);
            ++_deleted_leaf_number;
        }
        if(_nb_lig <= 1) {
            _kill_culm = true;
            if(_index == 1) {
            } else {
                _deleted_realloc_biomass = (_realocationCoeff) * (_internode_biomass_sum + _leaf_biomass_sum);
                _deleted_senesc_dw = (1 - _realocationCoeff) * (_internode_biomass_sum + _leaf_biomass_sum);
            }

            //kill leaves and internodes
            std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
            while(it != _phytomer_models.end()) {
                (*it)->internode()->culm_dead(t);
                (*it)->kill_leaf(t);
                ++it;
            }

            //realocation to plant_supply
            _realloc_supply = _deleted_realloc_biomass + _culm_deficit + _culm_stock;


        }
    }

    double get_leaf_biomass(double t, int index) const
    {
        double biomass = _phytomer_models[index]->leaf()->get < double >(t, LeafModel::BIOMASS);
        return biomass;
    }

    double get_leaf_blade_area(double t, int index) const
    {
        double blade_area = _phytomer_models[index]->leaf()->get < double >(t, LeafModel::BLADE_AREA);
        return blade_area;
    }

    int  get_first_ligulated_leaf_index(double t) const
    {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        int i = 0;
        while (it != _phytomer_models.end()) {
            if (not (*it)->is_leaf_dead(t) and (*it)->is_leaf_lig(t)) {
                break;
            }
            ++it;
            ++i;
        }
        if (it != _phytomer_models.end()) {
            return i;
        } else {
            return -1;
        }
    }

    int  get_first_alive_leaf_index(double t) const
    {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        int i = 0;
        int index = -1;
        while (it != _phytomer_models.end()) {
            if (not (*it)->is_leaf_dead(t) and ((*it)->leaf()->get < double >(t, LeafModel::BIOMASS) > 0)) {
                index = i;
                break;
            }
            ++it;
            ++i;
        }
        return index;
    }

    int  get_first_alive_leaf_index2(double t) const
    {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        int i = 0;
        int index = -1;
        while (it != _phytomer_models.end()) {
            if (not (*it)->is_leaf_dead(t-1) and ((*it)->leaf()->get < double >(t-1, LeafModel::BIOMASS) > 0)) {
                index = i;
                break;
            }
            ++it;
            ++i;
        }
        return index;
    }

    int  get_first_alive_leaf_creation_date(double t) const
    {
        std::deque < PhytomerModel* >::const_iterator it = _phytomer_models.begin();
        double creation_date = -1;
        int i = 1;
        while (it != _phytomer_models.end()) {
            if (not (*it)->is_leaf_dead(t) and ((*it)->leaf()->get < double >(t, LeafModel::BIOMASS) > 0)) {
                creation_date = (*it)->leaf()->get < double >(t, LeafModel::FIRST_DAY);
                break;
            }
            ++it;
            i++;
        }
        return creation_date;
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        _nb_leaf_param2 = _parameters.get("nb_leaf_param2");
        _plasto_init = _parameters.get("plasto_init");
        _phyllo_init = _parameters.get("phyllo_init");
        _ligulo_init = _parameters.get("ligulo_init");
        _coeff_Plasto_PI = _parameters.get("coef_plasto_PI");
        _coeff_Phyllo_PI = _parameters.get("coef_phyllo_PI");
        _coeff_Ligulo_PI = _parameters.get("coef_ligulo_PI");

        if(_plant_phenostage >= _nb_leaf_param2) {
            _plasto = _plasto_init * _coeff_Plasto_PI;
            _phyllo = _phyllo_init * _coeff_Phyllo_PI;
            _ligulo = _ligulo_init * _coeff_Ligulo_PI;
        } else {
            _plasto = _plasto_init;
            _phyllo = _phyllo_init;
            _ligulo = _ligulo_init;
        }
        PhytomerModel* first_phytomer = new PhytomerModel(1, _is_first_culm, _plasto, _phyllo, _ligulo, _LL_BL);
        setsubmodel(PHYTOMERS, first_phytomer);
        first_phytomer->init(t, parameters);
        _phytomer_models.push_back(first_phytomer);
        if(_is_first_culm) {
            PhytomerModel* second_phytomer = new PhytomerModel(2, _is_first_culm, _plasto, _phyllo, _ligulo, _LL_BL);
            PhytomerModel* third_phytomer = new PhytomerModel(3, _is_first_culm, _plasto, _phyllo, _ligulo, _LL_BL);
            PhytomerModel* fourth_phytomer = new PhytomerModel(4, _is_first_culm, _plasto, _phyllo, _ligulo, _LL_BL);

            setsubmodel(PHYTOMERS, second_phytomer);
            second_phytomer->init(t, parameters);
            _phytomer_models.push_back(second_phytomer);

            setsubmodel(PHYTOMERS, third_phytomer);
            third_phytomer->init(t, parameters);
            _phytomer_models.push_back(third_phytomer);

            setsubmodel(PHYTOMERS, fourth_phytomer);
            fourth_phytomer->init(t, parameters);
            _phytomer_models.push_back(fourth_phytomer);
        }
        _culm_stock_model->init(t, parameters);
        _culm_ictmodel->init(t, parameters);
        _culm_thermaltime_model->put(t, ThermalTimeModel::IS_FIRST_CULM, _is_first_culm);
        _culm_thermaltime_model->init(t, parameters);
        _culm_thermaltime_modelNG->init(t, parameters);

        last_time = t-1;
        first_day = t;

        //parameters
        _nb_leaf_max_after_pi = _parameters.get("nb_leaf_max_after_PI");
        _phenostage_pre_flo_to_flo  = _parameters.get("phenostage_PRE_FLO_to_FLO");
        _coeff_pi_lag = _parameters.get("coeff_PI_lag");
        _realocationCoeff = _parameters.get("realocationCoeff");
        _maxleaves = _parameters.get("maxleaves");

        //    internals
        _tt_plasto = _plasto_init;
        _tt_phyllo = _phyllo_init;
        _tt_ligulo = _ligulo_init;
        _nb_lig = 0;
        _stem_leaf_predim = 0;
        _leaf_biomass_sum = 0;
        _leaf_last_demand_sum = 0;
        _leaf_demand_sum = 0;
        _internode_last_demand_sum = 0;
        _internode_demand_sum = 0;
        _internode_biomass_sum = 0;
        _internode_len_sum = 0;
        _leaf_blade_area_sum = 0;
        _last_ligulated_leaf = -1;
        _last_ligulated_leaf_len = 0;
        _realloc_biomass_sum = 0;
        _senesc_dw_sum = 0;
        _last_leaf_biomass_sum = 0;
        _panicle_day_demand = 0;
        _panicle_weight = 0;
        _peduncle_inerlen_predim = 0;
        _peduncle_inerdiam_predim = 0;
        _peduncle_biomass = 0;
        _peduncle_last_demand = 0;
        _peduncle_day_demand = 0;
        _peduncle_len = 0;
        _first_leaf_len = 0;
        _deleted_leaf_number = 0;
        _deleted_senesc_dw = 0;
        _nb_lig_tot = 0;
        _kill_culm = false;
        _last_ligulated_leaf_blade_len = 0;
        _deleted_realloc_biomass = 0;
        _culm_test_ic = 0;
        _culm_stock = 0;
        _culm_deficit = 0;
        _last_leaf_index = -1;
        _realloc_supply = 0;
        _creation_date = t;
        _plant_biomass_sum = 0;
        _plant_blade_area_sum = 0;
        _plant_leaf_biomass_sum = 0;
        _last_plant_biomass_sum = 0;
        _started_PI = false;
        _culm_phase = culm::INITIAL;
        _last_phase = _culm_phase;
        _culm_phenostage = 1;
        _culm_appstage = 0;
        _culm_ligstage = 0;
        _culm_appstage_cste = _culm_appstage;
        _culm_ligstage_cste = _culm_ligstage;
        _culm_phenostage_cste = _culm_phenostage;
        _culm_phenostage_at_lag = 0;
        _culm_phenostage_at_pi = 0;
        _culm_phenostage_at_pre_flo = 0;
        _lag = false;
        _is_lagged = false;
        _culm_EDD = 0;
        _culm_DD = 0;
        _culm_bool_crossed_plasto = 0;
        _last_leaf_blade_area = 0;
        _nb_leaf = 1;
        _nb_init_leaves = 0;
        _nb_app_leaves = 0;
        _nb_app_leaves_tot = 0;
        _leaf_delay = 0;
        _reductionLER = 1.;
        _sheath_LLL = 0;
        _bool_crossed_phyllo = -1;
        _bool_crossed_ligulo = -1;
        _culm_bool_crossed_phyllo = -1;
        _culm_bool_crossed_ligulo = -1;
        _lig_index = 0;
        _culm_nbleaf_param_2 = _nb_leaf_param2;
        _culm_plasto_visu = _plasto_init;
        _culm_ligulo_visu = _ligulo_init;
        _culm_phyllo_visu = _phyllo_init;
        _culm_DD_phyllo = 0;
        _culm_DD_ligulo = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;

    //  submodels
    std::unique_ptr < CulmStockModelNG > _culm_stock_model;
    std::unique_ptr < IctModel > _culm_ictmodel;
    std::unique_ptr < ThermalTimeModel > _culm_thermaltime_model;
    std::unique_ptr < ThermalTimeModelNG > _culm_thermaltime_modelNG;
    std::deque < PhytomerModel* > _phytomer_models;
    std::unique_ptr < PanicleModel > _panicle_model;
    std::unique_ptr < PeduncleModel > _peduncle_model;

    //    attributes
    double _index;
    bool _is_first_culm;


    //parameters
    double _nb_leaf_max_after_pi;
    double _phenostage_pre_flo_to_flo;
    double _coeff_pi_lag;
    double _nb_leaf_param2;
    double _realocationCoeff;
    double _maxleaves;
    double _plasto_init;
    double _ligulo_init;
    double _phyllo_init;
    double _coeff_Plasto_PI;
    double _coeff_Phyllo_PI;
    double _coeff_Ligulo_PI;

    double first_day;

    //    internals
    double _phyllo;
    double _ligulo;
    double _tt_plasto;
    double _tt_phyllo;
    double _tt_ligulo;
    bool _started_PI;
    double _culm_phenostage;
    double _culm_appstage;
    double _culm_ligstage;
    double _culm_appstage_cste;
    double _culm_ligstage_cste;
    double _culm_phenostage_at_lag;
    culm::culm_phase _culm_phase;
    culm::culm_phase _last_phase;
    bool _lag;
    bool _is_lagged;
    double _nb_lig;
    double _stem_leaf_predim;
    double _leaf_biomass_sum;
    double _leaf_last_demand_sum;
    double _leaf_demand_sum;
    double _internode_last_demand_sum;
    double _internode_demand_sum;
    double _internode_biomass_sum;
    double _internode_len_sum;
    double _leaf_blade_area_sum;
    double _realloc_biomass_sum;
    double _senesc_dw_sum;
    int _last_ligulated_leaf;
    double _last_ligulated_leaf_len;
    double _last_leaf_biomass_sum;
    double _panicle_day_demand;
    double _panicle_weight;
    double _peduncle_inerlen_predim;
    double _peduncle_inerdiam_predim;
    double _peduncle_biomass;
    double _peduncle_last_demand;
    double _peduncle_day_demand;
    double _peduncle_len;
    double _culm_phenostage_at_pi;
    double _culm_phenostage_at_pre_flo;
    double _first_leaf_len;
    double _deleted_leaf_number;
    double _deleted_senesc_dw;
    double _nb_lig_tot;
    bool _kill_culm;
    double _last_ligulated_leaf_blade_len;
    double _deleted_realloc_biomass;
    double _culm_test_ic;
    double _culm_stock;
    double _culm_deficit;
    double _last_leaf_index;
    double _realloc_supply;
    int _culm_phenostage_cste;
    double _culm_EDD;
    double _culm_DD;
    double _culm_bool_crossed_plasto;
    double _creation_date;
    double _last_leaf_blade_area;
    double _nb_leaf;
    double _leaf_delay;
    double _reductionLER;
    double _sheath_LLL;
    double _nb_init_leaves;
    double _nb_app_leaves;
    double _nb_app_leaves_tot;
    double _culm_bool_crossed_phyllo;
    double _culm_bool_crossed_ligulo;
    double _lig_index;
    double _culm_nbleaf_param_2;
    double _culm_plasto_visu;
    double _culm_ligulo_visu;
    double _culm_phyllo_visu;
    double _culm_DD_phyllo;
    double _culm_DD_ligulo;

    //    externals
    int _plant_phenostage;
    int _plant_appstage;
    int _plant_ligstage;
    plant::plant_state _plant_state;
    plant::plant_phase _plant_phase;
    double _MGR;
    double _LL_BL;
    double _plasto;
    double _bool_crossed_plasto;
    double _bool_crossed_phyllo;
    double _bool_crossed_ligulo;
    double _dd;
    double _edd;
    double _delta_t;
    double _ftsw;
    double _fcstr;
    double _predim_leaf_on_mainstem;
    double _sla;
    double _test_ic;
    double _plant_stock;
    double _plant_deficit;
    double _plant_biomass_sum;
    double _plant_leaf_biomass_sum;
    double _plant_blade_area_sum;
    double _assim;
    double _last_plant_biomass_sum;
    bool _is_first_day_pi;
    double _ms_sheath_LLL;
    int _ms_phyt_index;
};

} // namespace model
#endif
