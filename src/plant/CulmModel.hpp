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

#include <defines.hpp>
#include <plant/processes/CulmStockModel.hpp>
#include <plant/processes/ThermalTimeModelNG.hpp>
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
                     CULM_DEFICIT, CULM_STOCK, LAST_LEAF_INDEX, REALLOC_SUPPLY, PANICLE_WEIGHT, LAST_LEAF_BLADE_AREA, NBLEAF };

    enum externals { BOOL_CROSSED_PLASTO, DD, EDD, DELTA_T, FTSW, FCSTR, PHENO_STAGE,
                     PREDIM_LEAF_ON_MAINSTEM, SLA, PLANT_PHASE,
                     PLANT_STATE, TEST_IC, PLANT_STOCK,
                     PLANT_DEFICIT, ASSIM, MGR, PLASTO, LIGULO, LL_BL,IS_FIRST_DAY_PI,
                   };


    CulmModel(int index):
        _index(index), _is_first_culm(index == 1),
        _culm_stock_model(new CulmStockModelNG),
        _culm_thermaltime_model(new ThermalTimeModelNG)
    {
        // submodels
        //        Submodels( ((CULM_STOCK, _culm_stock_model.get())) );
        //        @TODO gérer le problème de sous-modèle non calculé pour culmstock

        //    internals

        //@TODO regarder pourquoi ca pointe dans le vide
        //        InternalS(STOCK, _culm_stock_model.get(), CulmStockModel::STOCK);
        //        InternalS(DEFICIT, _culm_stock_model.get(), CulmStockModel::DEFICIT);
        //        InternalS(SURPLUS, _culm_stock_model.get(), CulmStockModel::SURPLUS);


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

        //    externals
        External(BOOL_CROSSED_PLASTO, &CulmModel::_bool_crossed_plasto);
        External(DD, &CulmModel::_dd);
        External(EDD, &CulmModel::_edd);
        External(DELTA_T, &CulmModel::_delta_t);
        External(FTSW, &CulmModel::_ftsw);
        External(FCSTR, &CulmModel::_fcstr);
        External(PHENO_STAGE, &CulmModel::_plant_phenostage);
        External(PREDIM_LEAF_ON_MAINSTEM, &CulmModel::_predim_leaf_on_mainstem);
        External(SLA, &CulmModel::_sla);
        External(PLANT_STATE, &CulmModel::_plant_state);
        External(PLANT_PHASE, &CulmModel::_plant_phase);
        External(TEST_IC, &CulmModel::_test_ic);
        External(PLANT_STOCK, &CulmModel::_plant_stock);
        External(PLANT_DEFICIT, &CulmModel::_plant_deficit);
        External(ASSIM, &CulmModel::_assim);
        External(MGR, &CulmModel::_MGR);
        External(PLASTO, &CulmModel::_plasto);
        External(LIGULO, &CulmModel::_ligulo);
        External(LL_BL, &CulmModel::_LL_BL);
        External(IS_FIRST_DAY_PI, &CulmModel::_is_first_day_pi);


    }

    virtual ~CulmModel()
    {
        _culm_stock_model.reset(nullptr);
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
                && (get_phytomer_number() < _nb_leaf_pi + _nb_leaf_max_after_pi)
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
            if (_culm_phenostage == _culm_phenostage_at_pi + _nb_leaf_max_after_pi + 1) {
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
            if(_culm_phenostage == _culm_phenostage_at_pre_flo + _phenostage_pre_flo_to_flo ) {
                _culm_phase = culm::FLO;
            }
            break;
        }
        }
    }

    void compute(double t, bool /* update */) {
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
            return;
        }

        //thermaltimevalues update
        //qDebug() << "test0";
        if(_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE or _is_first_day_pi or _creation_date == t) {
            _culm_EDD = _edd;
            _culm_DD = _dd;
            _culm_bool_crossed_plasto = _bool_crossed_plasto;
            if(_culm_bool_crossed_plasto >= 0) {
                _culm_phenostage = _culm_phenostage +1;
                _culm_phenostage_cste = _culm_phenostage;
            }
        } else {
            //qDebug() << "test1";
            _culm_DD = _culm_thermaltime_model->get < double >(t, ThermalTimeModelNG::CULM_DD);
            _culm_EDD = _culm_thermaltime_model->get < double >(t, ThermalTimeModelNG::EDD);
            _culm_bool_crossed_plasto = _culm_thermaltime_model->get < double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PLASTO);
            _culm_phenostage = _culm_phenostage_cste + _culm_thermaltime_model->get < int >(t, ThermalTimeModelNG::PHENO_STAGE);
            //qDebug() << "test2";
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

        //Nécessaire pour coder l'erreur delphi du passage de l'entrenoeud à réalisation à partir d'un certain index
        if(_culm_phase == culm::VEGETATIVE or _culm_phase == culm::INITIAL) {
            _last_leaf_index = _phytomer_models.size();
        }

        step_state(t);

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


        while (it != _phytomer_models.end()) {
            //Phytomers
            compute_phytomers(it, previous_it, i, t);
            //Sum
            compute_vars(it, previous_it, i, t);
            //GetLastINnonVegetative
            get_nonvegetative_in(it, t);

            if(i == 0 and not (*it)->is_leaf_dead()) {
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
            (*_peduncle_model)(t);
        }

        if(_peduncle_model.get()) {
            _peduncle_last_demand = _peduncle_model->get < double >(t, PeduncleModel::LAST_DEMAND);
            _peduncle_day_demand = _peduncle_model->get < double >(t, PeduncleModel::DEMAND);
            _peduncle_len = _peduncle_model->get < double >(t, PeduncleModel::LENGTH);
            _peduncle_biomass = _peduncle_model->get < double >(t, PeduncleModel::BIOMASS);
        }

        // @TODO : Ajouter un test à partir de ce moment pour tuer les talles sans feuilles du à senesc
        //        if(_culm_phase != culm::INITIAL and _culm_phase != culm::VEGETATIVE) {
        //            delete_culmNG(t);
        //        }
    }

    void get_nonvegetative_in(std::deque < PhytomerModel* >::iterator it, double t) {
        if(_peduncle_model.get()) {
            //@TODO : à modifier pour prendre la phase de l'entrenoeud ?
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

    void compute_thermaltime(double t) {
        _culm_thermaltime_model->put(t, ThermalTimeModelNG::CULM_STOCK, _culm_stock);
        _culm_thermaltime_model->put(t, ThermalTimeModelNG::CULM_DEFICIT, _culm_deficit);
        _culm_thermaltime_model->put(t, ThermalTimeModelNG::PLASTO, _plasto);
        _culm_thermaltime_model->put(t, ThermalTimeModelNG::IS_FIRST_DAY_OF_INDIVIDUALIZATION, _is_first_day_pi);
        _culm_thermaltime_model->put(t, ThermalTimeModelNG::PLASTO_DELAY, 0.); //@TODO : mettre vraie valeur de plasto_delay
        (*_culm_thermaltime_model)(t);
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
        (*it)->internode()->put(t, InternodeModel::CULM_DEFICIT, _culm_deficit);
        (*it)->internode()->put(t, InternodeModel::CULM_STOCK, _culm_stock);
        (*it)->internode()->put(t, InternodeModel::CULM_PHASE, _culm_phase);
        (*it)->internode()->put(t, InternodeModel::PLASTO, _plasto);
        (*it)->internode()->put(t, InternodeModel::LIGULO, _ligulo);
        (*it)->internode()->put(t, InternodeModel::NB_LIG, _nb_lig);
        (*it)->internode()->put(t, InternodeModel::BOOL_CROSSED_PLASTO, _culm_bool_crossed_plasto);
        (*it)->internode()->put(t, InternodeModel::LAST_LEAF_INDEX, _last_leaf_index);
        (*it)->internode()->put(t, InternodeModel::PHENOSTAGE, _plant_phenostage);
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
        if (not (*it)->is_leaf_dead()) {
            if((*it)->is_leaf_lig(t)) {
                ++ _nb_lig;
                ++ _nb_lig_tot;
                double ratio = (*it)->leaf()->get < double >(t, LeafModel::BLADE_AREA) / (*it)->leaf()->get < double >(t, LeafModel::LAST_BLADE_AREA);
                if(ratio >= 0.5) {
                    _nb_leaf = _nb_leaf +1;
                }
            } else {
                _nb_leaf = _nb_leaf + 1;
            }

            if (_index == 1) {
                _stem_leaf_predim = (*it)->get < double, LeafModel >(t, PhytomerModel::LEAF_PREDIM);
                if((*it)->is_leaf_lig(t) and not( (*it)->is_leaf_dead())) {
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

            PhytomerModel* phytomer = new PhytomerModel(index, _is_first_culm, _plasto, _ligulo, _LL_BL);
            setsubmodel(PHYTOMERS, phytomer);
            phytomer->init(t, _parameters);
            _phytomer_models.push_back(phytomer);
        }
    }

    int get_phytomer_number() const
    { return _phytomer_models.size(); }

    int get_alive_phytomer_number() const
    { return _phytomer_models.size() - _deleted_leaf_number; }

    CulmStockModelNG * stock_model() const
    { return _culm_stock_model.get(); }

    ThermalTimeModelNG * thermaltime_model() const
    { return _culm_thermaltime_model.get(); }

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
        std::string date = artis::utils::DateTime::toJulianDayFmt(t, artis::utils::DATE_FORMAT_YMD);
        if(index != -1) {
            _deleted_senesc_dw = (1 - _realocationCoeff) * _phytomer_models[index]->leaf()->get < double >(t-1, LeafModel::BIOMASS);
            _deleted_realloc_biomass = (_realocationCoeff) * _phytomer_models[index]->leaf()->get < double >(t-1, LeafModel::BIOMASS);
            _phytomer_models[index]->kill_leaf(t);
            //qDebug() << "Le : " << QString::fromStdString(date) << " on tue la feuille" << index + 1 << " sur la talle " << _index - 1;
            ++_deleted_leaf_number;
        }
        if(_nb_lig <= 1) {
            _kill_culm = true;
            if(_index == 1) {
                //qDebug() << "Plus de feuilles lig sur le brin maître, la plante est morte";
            } else {
                _deleted_realloc_biomass = (_realocationCoeff) * (_internode_biomass_sum + _leaf_biomass_sum);
                _deleted_senesc_dw = (1 - _realocationCoeff) * (_internode_biomass_sum + _leaf_biomass_sum);
                //qDebug() << "Plus de feuilles lig sur la talle " << _index - 1 << " on tue la talle";
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
            if (not (*it)->is_leaf_dead() and (*it)->is_leaf_lig(t)) {
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
            if (not (*it)->is_leaf_dead() and ((*it)->leaf()->get < double >(t, LeafModel::BIOMASS) > 0)) {
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
            if (not (*it)->is_leaf_dead() and ((*it)->leaf()->get < double >(t-1, LeafModel::BIOMASS) > 0)) {
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
            if (not (*it)->is_leaf_dead() and ((*it)->leaf()->get < double >(t, LeafModel::BIOMASS) > 0)) {
                creation_date = (*it)->leaf()->get < double >(t, LeafModel::FIRST_DAY);
                break;
            }
            ++it;
            i++;
        }
        return creation_date;
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        PhytomerModel* first_phytomer = new PhytomerModel(1, _is_first_culm, _plasto, _ligulo, _LL_BL);

        setsubmodel(PHYTOMERS, first_phytomer);
        first_phytomer->init(t, parameters);
        _phytomer_models.push_back(first_phytomer);
        _culm_stock_model->init(t, parameters);
        _culm_thermaltime_model->init(t, parameters);


        last_time = t-1;
        //parameters
        _parameters = parameters;
        _nb_leaf_pi = _parameters.get("nbleaf_pi");
        _nb_leaf_max_after_pi = _parameters.get("nb_leaf_max_after_PI");
        _phenostage_pre_flo_to_flo  = _parameters.get("phenostage_PRE_FLO_to_FLO");
        _coeff_pi_lag = _parameters.get("coeff_PI_lag");
        _nb_leaf_param2 = _parameters.get("nb_leaf_param2");
        _realocationCoeff = _parameters.get("realocationCoeff");

        //    internals
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

        _started_PI = false;
        _culm_phase = culm::INITIAL;
        _last_phase = _culm_phase;
        _culm_phenostage = 1;
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
    }

private:
    ecomeristem::ModelParameters _parameters;

    //  submodels
    std::unique_ptr < CulmStockModelNG > _culm_stock_model;
    std::unique_ptr < ThermalTimeModelNG > _culm_thermaltime_model;
    std::deque < PhytomerModel* > _phytomer_models;
    std::unique_ptr < PanicleModel > _panicle_model;
    std::unique_ptr < PeduncleModel > _peduncle_model;

    //    attributes
    double _index;
    bool _is_first_culm;


    //parameters
    double _nb_leaf_pi;
    double _nb_leaf_max_after_pi;
    double _phenostage_pre_flo_to_flo;
    double _coeff_pi_lag;
    double _nb_leaf_param2;
    double _realocationCoeff;

    //    internals
    bool _started_PI;
    double _culm_phenostage;
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

    //    externals
    int _plant_phenostage;
    plant::plant_state _plant_state;
    plant::plant_phase _plant_phase;
    double _MGR;
    double _LL_BL;
    double _plasto;
    double _ligulo;
    double _bool_crossed_plasto;
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
};

} // namespace model
#endif
