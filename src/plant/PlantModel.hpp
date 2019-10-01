/**
 * @file ecomeristem/plant/Model.hpp
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

#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <defines.hpp>
#include <plant/CulmModel.hpp>
#include <plant/RootModel.hpp>
#include <plant/processes/WaterBalanceModel.hpp>
#include <plant/processes/PlantStockModel.hpp>
#include <plant/processes/AssimilationModel.hpp>
#include <plant/processes/InterceptionModel.hpp>

using namespace model;

class PlantModel : public CoupledModel < PlantModel >
{
public:
    enum submodels { WATER_BALANCE, STOCK, ASSIMILATION, INTERCEPTION,
                     ROOT, CULMS};

    enum internals { LIG, APP, LEAF_BIOMASS_SUM, INTERNODE_BIOMASS_SUM,
                     SENESC_DW_SUM, LEAF_LAST_DEMAND_SUM,
                     INTERNODE_LAST_DEMAND_SUM, LEAF_DEMAND_SUM,
                     INTERNODE_DEMAND_SUM, PANICLE_DEMAND_SUM,
                     PLANT_PHASE, PLANT_STATE, PAI, HEIGHT, HEIGHT_P,
                     PLASTO_INIT, PHYLLO_INIT, LIGULO_INIT, TT_LIG, IH,
                     BIOMLEAF, BIOMIN, BIOMINSTRUCT, INTERNODE_STOCK_SUM,
                     REALLOC_BIOMASS_SUM, PEDUNCLE_BIOMASS_SUM, PEDUNCLE_LAST_DEMAND_SUM,
                     CULM_SURPLUS_SUM, QTY, LL_BL, PLANT_STOCK, REALLOC_SUM_SUPPLY,
                     TA, DELTA_T, TT, BOOL_CROSSED_PLASTO, BOOL_CROSSED_PHYLLO, BOOL_CROSSED_LIGULO,
                     EDD, DD, PHENOSTAGE, APPSTAGE, LIGSTAGE, SLA, DD_PHYLLO, DD_LIGULO, PHYLLO_VISU,
                     TILLERNB_1, NBLEAF, BIOMAERO2, BIOMLEAFMAINSTEM, BIOMINMAINSTEM, AREALFEL,
                     MAINSTEM_STOCK_IN, BIOMINMAINSTEMSTRUCT, BIOMLEAFMAINSTEMSTRUCT, MAINSTEM_STOCK,
                     DEAD_LEAF_NB, INTERNODE_LENGTH_MAINSTEM, PANICLE_MAINSTEM_DW,
                     PANICLE_DW, LEAF_DELAY, PHENOSTAGE_AT_FLO, LIG_INDEX,
                     MS_INDEX, DELETED_LEAF_BIOMASS, VISI, PREDIM_APP_LEAF_MS, NB_CR_TILLERS, TAE, PANICLENB,
                     TOTAL_LENGTH_MAINSTEM, BIOMMAINSTEM, SLAPLANT, BIOMLEAFTOT, SENESC_DW, BIOMINSHEATHMS,
                     BIOMINSHEATH, BIOMAEROTOT, INTERC1, INTERC2, MS_LEAF2_LEN, PARI, BIOMAEROFW, TILLERFW,
                     MAINSTEMFW, MAINSTEMBLADEFW, TILLERLEAFFW, FIRST_DAY_INDIV,
                     CREATED_TILLERS, CULM_DEFICIT_SUM/*, CURRENT_DATE*/ };

    PlantModel() :
        _water_balance_model(new WaterBalanceModel),
        _stock_model(new PlantStockModel),
        _assimilation_model(new AssimilationModel),
        _interception_model(new InterceptionModel),
        _root_model(new RootModel)
    {
        // submodels
        subModel(WATER_BALANCE, _water_balance_model.get());
        subModel(STOCK, _stock_model.get());
        subModel(ASSIMILATION, _assimilation_model.get());
        subModel(INTERCEPTION, _interception_model.get());
        subModel(ROOT, _root_model.get());

        // local internals
        Internal( LIG, &PlantModel::_lig );
        Internal( APP, &PlantModel::_app );
        Internal( LEAF_BIOMASS_SUM, &PlantModel::_leaf_biomass_sum );
        Internal( LEAF_DEMAND_SUM, &PlantModel::_leaf_demand_sum );
        Internal( LEAF_LAST_DEMAND_SUM, &PlantModel::_leaf_last_demand_sum );
        Internal( INTERNODE_BIOMASS_SUM, &PlantModel::_internode_biomass_sum );
        Internal( INTERNODE_DEMAND_SUM, &PlantModel::_internode_demand_sum );
        Internal( INTERNODE_LAST_DEMAND_SUM, &PlantModel::_internode_last_demand_sum );
        Internal( PANICLE_DEMAND_SUM, &PlantModel::_panicle_demand_sum );
        Internal( SENESC_DW_SUM, &PlantModel::_senesc_dw_sum );
        Internal( PAI, &PlantModel::_leaf_blade_area_sum );
        Internal( HEIGHT, &PlantModel::_height );
        Internal( HEIGHT_P, &PlantModel::_height_ped );
        Internal( PLANT_PHASE, &PlantModel::_plant_phase );
        Internal( PLANT_STATE, &PlantModel::_plant_state );
        Internal( PLASTO_INIT, &PlantModel::_plasto_init );
        Internal( PHYLLO_INIT, &PlantModel::_phyllo_init );
        Internal( LIGULO_INIT, &PlantModel::_ligulo_init );
        Internal( TT_LIG, &PlantModel::_TT_lig );
        Internal( IH, &PlantModel::_IH );
        Internal( BIOMLEAF, &PlantModel::_biomLeaf );
        Internal( BIOMIN, &PlantModel::_biomin );
        Internal( BIOMINSTRUCT, &PlantModel::_biominstruct );
        Internal( INTERNODE_STOCK_SUM, &PlantModel::_internode_stock_sum );
        Internal( REALLOC_BIOMASS_SUM, &PlantModel::_realloc_biomass_sum );
        Internal( PEDUNCLE_BIOMASS_SUM, &PlantModel::_peduncle_biomass_sum );
        Internal( PEDUNCLE_LAST_DEMAND_SUM, &PlantModel::_peduncle_last_demand_sum );
        Internal( CULM_SURPLUS_SUM, &PlantModel::_culm_surplus_sum );
        Internal( QTY, &PlantModel::_qty );
        Internal( LL_BL, &PlantModel::_LL_BL );
        Internal( PLANT_STOCK, &PlantModel::_stock );
        Internal( REALLOC_SUM_SUPPLY, &PlantModel::_realloc_sum_supply );
        Internal( TA, &PlantModel::_Ta );
        Internal( DELTA_T, &PlantModel::_deltaT );
        Internal( TT, &PlantModel::_TT );
        Internal( BOOL_CROSSED_PLASTO, &PlantModel::_bool_crossed_plasto );
        Internal( BOOL_CROSSED_PHYLLO, &PlantModel::_bool_crossed_phyllo );
        Internal( BOOL_CROSSED_LIGULO, &PlantModel::_bool_crossed_ligulo );
        Internal( EDD, &PlantModel::_EDD );
        Internal( DD, &PlantModel::_DD );
        Internal( PHENOSTAGE, &PlantModel::_phenostage );
        Internal( APPSTAGE, &PlantModel::_appstage );
        Internal( LIGSTAGE, &PlantModel::_ligstage );
        Internal( SLA, &PlantModel::_sla );
        Internal( TILLERNB_1, &PlantModel::_tillerNb_1 );
        Internal( NBLEAF, &PlantModel::_nbleaf );
        Internal( BIOMAERO2, &PlantModel::_biomAero2 );
        Internal( BIOMLEAFMAINSTEM, &PlantModel::_biomLeafMainstem );
        Internal( BIOMINMAINSTEM, &PlantModel::_biomInMainstem );
        Internal( AREALFEL, &PlantModel::_areaLFEL );
        Internal( MAINSTEM_STOCK_IN, &PlantModel::_mainstem_stock_IN );
        Internal( BIOMINMAINSTEMSTRUCT, &PlantModel::_biomInMainstemstruct );
        Internal( BIOMLEAFMAINSTEMSTRUCT, &PlantModel::_biomLeafMainstemstruct );
        Internal( MAINSTEM_STOCK, &PlantModel::_mainstem_stock );
        Internal( DEAD_LEAF_NB, &PlantModel::_deadleafNb );
        Internal( INTERNODE_LENGTH_MAINSTEM, &PlantModel::_internode_length_mainstem );
        Internal( PANICLE_MAINSTEM_DW, &PlantModel::_panicleMainstemDW );
        Internal( PANICLE_DW, &PlantModel::_panicleDW );
        Internal( LEAF_DELAY, &PlantModel::_leaf_delay );
        Internal( DD_PHYLLO, &PlantModel::_DD_phyllo );
        Internal( DD_LIGULO, &PlantModel::_DD_ligulo );
        Internal( PHYLLO_VISU, &PlantModel::_phyllo_visu );
        Internal( PHENOSTAGE_AT_FLO, &PlantModel::_phenostage_at_flo);
        Internal( LIG_INDEX, &PlantModel::_lig_index);
        Internal( MS_INDEX, &PlantModel::_ms_index);
        Internal( DELETED_LEAF_BIOMASS, &PlantModel::_deleted_leaf_biomass);
        Internal( VISI, &PlantModel::_visi);
        Internal( PREDIM_APP_LEAF_MS, &PlantModel::_predim_app_leaf_on_mainstem);
        Internal( NB_CR_TILLERS, &PlantModel:: _nb_tillers);
        Internal( TAE, &PlantModel:: _tae);
        Internal( PANICLENB, &PlantModel:: _paniclenb);
        Internal( TOTAL_LENGTH_MAINSTEM, &PlantModel::_total_length_mainstem);
        Internal( BIOMMAINSTEM, &PlantModel::_biomMainstem);
        Internal( SLAPLANT, &PlantModel::_slaplant);
        Internal( BIOMLEAFTOT, &PlantModel::_biomLeafTot);
        Internal( SENESC_DW, &PlantModel::_senesc_dw);
        Internal( BIOMINSHEATHMS, &PlantModel::_biomInSheathMainstem);
        Internal( BIOMINSHEATH, &PlantModel::_biomInSheath);
        Internal( BIOMAEROTOT, &PlantModel::_biomAeroTot);
        Internal( INTERC1, &PlantModel::_interc1);
        Internal( INTERC2, &PlantModel::_interc2);
        Internal( MS_LEAF2_LEN, &PlantModel::_ms_leaf2_len);
        Internal( PARI, &PlantModel::_pari);
        Internal( BIOMAEROFW, &PlantModel::_biomAeroFW);
        Internal( TILLERFW, &PlantModel::_tillerFW);
        Internal( MAINSTEMFW, &PlantModel::_biomMainstemFW);
        Internal( MAINSTEMBLADEFW, &PlantModel::_biomBladeMainstemFW);
        Internal( TILLERLEAFFW, &PlantModel::_tillerleafFW);
        Internal( FIRST_DAY_INDIV, &PlantModel::_is_first_day_pi);
        Internal( CREATED_TILLERS, &PlantModel::_createdTillers);
        Internal( CULM_DEFICIT_SUM, &PlantModel::_culm_deficit_sum);
        //Internal(CURRENT_DATE, &PlantModel::_current_date);


    }

    virtual ~PlantModel()
    {
        _root_model.reset(nullptr);
        _water_balance_model.reset(nullptr);
        _stock_model.reset(nullptr);
        _assimilation_model.reset(nullptr);
        _interception_model.reset(nullptr);

        auto it = _culm_models.begin();
        while (it != _culm_models.end()) {
            delete *it;
            ++it;
        }
    }


    bool is_phytomer_creatable() {
        return (_plant_phase == plant::VEGETATIVE
                || _plant_phase == plant::ELONG
                || _plant_phase == plant::PI
                || _plant_phase == plant::PRE_FLO
                );
    }

    void step_state(double t) {
        double ic = _stock_model->get <double> (t-1, PlantStockModel::IC);
        double FTSW  = _water_balance_model->get<double> (t, WaterBalanceModel::FTSW);
        if(_is_first_day_pi) {
            _is_first_day_pi = false;
        }

        if (FTSW <= 0 or ic <= -1) {
            // plant dead
            _plant_state << plant::NOGROWTH;
            return;
        }

        _plant_state >> plant::NEW_PHYTOMER_AVAILABLE;
        if (_stock <= 0) {
            _plant_state << plant::NOGROWTH;
        } else {
            _plant_state >> plant::NOGROWTH;
        }
        //Globals states
        if (_ligstage >= _nb_leaf_indiv) {
            if(!(_plant_state & plant::INDIV)) {
                _plant_state << plant::INDIV;
            }
            if(_ligstage == _nb_leaf_indiv & !_is_first_day_passed) {
                _is_first_day_pi = true;
                _is_first_day_passed = true;
            }
        }
        if(_stock <= 0) {
            return;
        }

        switch (_plant_phase) {
        case plant::INITIAL: {
            _plant_phase = plant::VEGETATIVE;
            break;
        }
        case plant::VEGETATIVE: {
            if(_ligstage == _nb_leaf_stem_elong and _phenostage < _maxleaves + 1) {
                _plant_phase = plant::ELONG;
                //_is_first_day_pi = true;
            } else if (_phenostage == _maxleaves + 1) {
                _plant_phase = plant::PI;
                //_is_first_day_pi = true;
            }
            break;
        }
        case plant::ELONG: {
            _is_first_day_pi = false;
            if(_phenostage == _maxleaves + 1) {
                _plant_phase = plant::PI;
            }
        }
        case plant::PI: {
            _is_first_day_pi = false;
            if (_ligstage == _maxleaves) {
                _plant_phase = plant::PRE_FLO;
            }
            break;
        }
        case plant::PRE_FLO: {
            if (_ligstage == _maxleaves + _phenostage_pre_flo_to_flo) {
                _plant_phase = plant::FLO;
                _phenostage_at_flo = _phenostage;
            }
            break;
        }
        case plant::FLO: {
            if (_phenostage == _phenostage_at_flo + _phenostage_to_end_filling) {
                _plant_phase = plant::END_FILLING;
            }
            break;
        }
        case plant::END_FILLING: {
            if (_phenostage == _phenostage_at_flo + _phenostage_to_end_filling + _phenostage_to_maturity) {
                _plant_phase = plant::MATURITY;
            }
            break;
        }}

        if ( _bool_crossed_plasto >= 0 && is_phytomer_creatable()) {
            _plant_state << plant::NEW_PHYTOMER_AVAILABLE;
        }
    }


    void compute(double t, bool /* update */) {
        //_current_date = t - _parameters.beginDate;
        //Delete culm and leaf
        std::deque < CulmModel* >::const_iterator nc = _culm_models.begin();
        while(nc != _culm_models.end()) {
            if((*nc)->get < bool, CulmModel >(t-1, CulmModel::KILL_CULM)) {
                _deleted_leaf_biomass += (*nc)->get < double, CulmModel >(t-1, CulmModel::DEL_LEAF_BIOM);
            }
            ++nc;
        }
        delete_leaf(t);

        if (_deleted_leaf_biomass > 0) {
            _qty = _deleted_leaf_biomass * _realocationCoeff;
            _stock = std::max(0., _qty + _stock_model->get < double >(t-1, PlantStockModel::DEFICIT));
            _deficit = std::min(0., _qty + _stock_model->get < double >(t-1, PlantStockModel::DEFICIT));
        } else {
            _stock = _stock_model->get < double >(t-1, PlantStockModel::STOCK);
            _deficit = _stock_model->get < double >(t-1, PlantStockModel::DEFICIT);
        }

        //Compute IC
        _stock_model->compute_IC(t);

        //Thermal time
        _Ta = _parameters.get(t).Temperature;
        _deltaT = _Ta - _Tb;
        //Thermal time New : TODO
        //if(_Ta < _Tb) {
        //    _deltaT = 0;
        //} else if(_Ta < _Topt) {
        //    _deltaT = _Ta - Tb;
        //} else if(_Ta < _Tmax) {
        //    _deltaT = _Tmax - _Ta;
        //} else {
        //    _deltaT = 0;
        //}

        _TT = _TT + _deltaT;


        std::deque < CulmModel* >::const_iterator culms = _culm_models.begin();
        int i = 0;
        while(culms != _culm_models.end()) {
            (*culms)->ictmodel()->put < double >(t, IctModel::IC_plant, _stock_model->get <double> (t-1, PlantStockModel::IC));
            (*culms)->compute_ictmodel(t);
            if(/*(_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE)*/ !(_plant_state & plant::INDIV)) {
                if(!(*culms)->get < bool, CulmModel >(t-1, CulmModel::KILL_CULM)) {
                    (*culms)->thermaltime_model()->put < double >(t, ThermalTimeModel::DELTA_T, _deltaT);
                    (*culms)->thermaltime_model()->put < double >(t, ThermalTimeModel::PLASTO_DELAY, _leaf_delay);
                    (*culms)->thermaltime_model()->put < double >(t, ThermalTimeModel::STOCK, _stock);
                    (*culms)->thermaltime_model()->put < plant::plant_state >(t, ThermalTimeModel::PLANT_STATE, _plant_state);
                    (*culms)->compute_thermaltime(t);
                    if(i == 0) {
                        _bool_crossed_plasto = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::BOOL_CROSSED_PLASTO);
                        _bool_crossed_phyllo = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::BOOL_CROSSED_PHYLLO);
                        _bool_crossed_ligulo = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::BOOL_CROSSED_LIGULO);
                        _EDD = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::EDD);
                        _DD = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::DD);
                        _ligulo_visu = (*culms)->thermaltime_model()->get< double >(t, ThermalTimeModel::LIGULO_VISU);
                        _phenostage = (*culms)->thermaltime_model()->get< int >(t, ThermalTimeModel::PHENO_STAGE);
                        _appstage = (*culms)->thermaltime_model()->get< int >(t, ThermalTimeModel::APP_STAGE);
                        _ligstage = (*culms)->thermaltime_model()->get< int >(t, ThermalTimeModel::LIG_STAGE);
                        _phenostage_cste = _phenostage;
                        _appstage_cste = _appstage;
                        _ligstage_cste = _ligstage;
                    }
                }
            } else {
                if(!(*culms)->get < bool, CulmModel >(t-1, CulmModel::KILL_CULM)) {
                    (*culms)->thermaltime_modelNG()->put < double >(t, ThermalTimeModelNG::DELTA_T, _deltaT);
                    (*culms)->thermaltime_modelNG()->put < double >(t, ThermalTimeModelNG::DELTA_T, _deltaT);
                    (*culms)->thermaltime_modelNG()->put < plant::plant_state >(t, ThermalTimeModelNG::PLANT_STATE, _plant_state);
                    (*culms)->compute_thermaltimeNG(t);
                    if(i == 0) {
                        _bool_crossed_plasto = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PLASTO);
                        _bool_crossed_phyllo = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PHYLLO);
                        _bool_crossed_ligulo = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_LIGULO);
                        _EDD = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::EDD);
                        _DD = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_DD);
                        _ligulo_visu = (*culms)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_LIGULO_VISU);
                        _phenostage = _phenostage_cste + (*culms)->thermaltime_modelNG()->get< int >(t, ThermalTimeModelNG::PHENO_STAGE);
                        _appstage = _appstage_cste + (*culms)->thermaltime_modelNG()->get< int >(t, ThermalTimeModelNG::APP_STAGE);
                        _ligstage = _ligstage_cste + (*culms)->thermaltime_modelNG()->get< int >(t, ThermalTimeModelNG::LIG_STAGE);
                    }
                }
            }
            ++culms;
            ++i;
        }

        //SLA
        _sla = _FSLA - _SLAp * std::log(_phenostage);

        //Water balance
        _water_balance_model->put < double >(t, WaterBalanceModel::INTERC,
                                             _assimilation_model->get < double >(t-1, AssimilationModel::INTERC));
        (*_water_balance_model)(t);

        // Manager
        //std::cout << t - _parameters.beginDate << std::endl;
        //std::cout << "Plant state :" << _plant_state << std::endl;
        //std::cout << "Plant phase :" << _plant_phase << std::endl;
        step_state(t);
        //std::cout << "Plant state :" << _plant_state << std::endl;
        //std::cout << "Plant phase :" << _plant_phase << std::endl;
        //if(_plant_phase == plant::MATURITY) {
            //std::cout << "FIRST DAY OF MATURITY : " << t - _parameters.beginDate << std::endl;
        //}

        //LLBL - MGR
        if (_phenostage == _nb_leaf_param2 and _bool_crossed_plasto >= 0 and _stock > 0) {
            _LL_BL = _LL_BL_init + _slope_LL_BL_at_PI;
            _MGR = _MGR * _coeff_MGR_PI;
        } else if (_phenostage > _nb_leaf_param2 and _bool_crossed_plasto > 0 and _phenostage <= _maxleaves) {
            _LL_BL = _LL_BL_init + _slope_LL_BL_at_PI * (_phenostage+1-_nb_leaf_param2);
        }

        //Tillering
        double ic = _stock_model->get < double >(t-1, PlantStockModel::IC);

        std::deque < CulmModel* >::const_iterator it = _culm_models.begin();
        _tae = 0;
        while(it != _culm_models.end()) {
            if(!(*it)->get < bool, CulmModel >(t-1, CulmModel::KILL_CULM) and (*it)->get < bool, CulmModel >(t-1, CulmModel::IS_COMPUTED)) {
                double bcphyllo = -1;
                if(/*_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE*/ !(_plant_state & plant::INDIV) or _is_first_day_pi) {
                    bcphyllo = (*it)->thermaltime_model()->get< double >(t, ThermalTimeModel::BOOL_CROSSED_PHYLLO);
                } else {
                    bcphyllo = (*it)->thermaltime_modelNG()->get< double >(t, ThermalTimeModelNG::CULM_BOOL_CROSSED_PHYLLO);
                }
                int potential_leaf = bcphyllo >= 0 ? 1 : 0;
                if ((*it)->get_app_phytomer_number(t) + potential_leaf >= _nbleaf_enabling_tillering) {
                    ++_tae;
                }
            }
            it++;
        }
        if (ic > _Ict) {
            _nb_tillers = _nb_tillers + _nbExistingTillers;
        }
        if (_bool_crossed_plasto > 0 and _nb_tillers >= 1) {
            _nb_tillers = std::min(_nb_tillers, _tae);
            if(_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE) {
                create_culm(t, _nb_tillers);
                _nbExistingTillers = _nbExistingTillers + _nb_tillers;
            } else if(_plant_phase == plant::ELONG or _plant_phase == plant::PI) {
                if(_culm_deficit_sum >= 0 /*and _createdTillers < 30*/) {
                    create_culm(t, _nb_tillers);
                    _nbExistingTillers = _nbExistingTillers + _nb_tillers;
                }
            }
        }

        //CulmModel

        compute_culms(t);

        //Lig update
        _lig_1 = _lig;
        std::deque < CulmModel* >::const_iterator mainstem = _culm_models.begin();
        _lig = (*mainstem)->get <double, CulmModel>(t, CulmModel::NB_LIG_TOT);
        _app = (*mainstem)->get < double, CulmModel >(t, CulmModel::NB_APP_LEAVES_TOT);
        _visi = (*mainstem)->get < double, CulmModel >(t, CulmModel::NB_APP_LEAVES);

        //TT_Lig
        if (t != _parameters.beginDate) {
            if (_lig_1 == _lig) {
                if (!(_plant_state & plant::NOGROWTH)) {
                    _TT_lig = _TT_lig + _EDD;
                }
            } else {
                _TT_lig = 0;
            }
        }

        //IH
        if (!(_plant_state & plant::NOGROWTH)) {
            _IH = _lig + std::min(1., _TT_lig / _ligulo_visu);
        }

        //Interception TESTS
        if(_intercmodel != 1) {
            _interception_model->put < double >(t, InterceptionModel::PAI, _leaf_blade_area_sum);
            (*_interception_model)(t);
            _assimilation_model->put < double >(t, AssimilationModel::EXT_INTERC,
                                                _interception_model->get < double >(t, InterceptionModel::INTERC)/_parameters.get(t).Par);
        } else {
            _assimilation_model->put < double >(t, AssimilationModel::EXT_INTERC,0);
        }

        //Assimilation
        _assimilation_model->put < double >(t, AssimilationModel::CSTR,
                                            _water_balance_model->get < double >(t, WaterBalanceModel::CSTR));
        _assimilation_model->put < double >(t, AssimilationModel::FCSTR,
                                            _water_balance_model->get < double >(t, WaterBalanceModel::FCSTR));
        _assimilation_model->put < double >(t, AssimilationModel::FCSTRA,
                                            _water_balance_model->get < double >(t, WaterBalanceModel::FCSTRA));
        _assimilation_model->put < double >(t, AssimilationModel::PAI, _leaf_blade_area_sum);
        _assimilation_model->put < double >(t, AssimilationModel::LEAFBIOMASS, _leaf_biomass_sum);
        _assimilation_model->put < double >(t, AssimilationModel::INTERNODEBIOMASS, _internode_biomass_sum);
        (*_assimilation_model)(t);

        _interc1 = _assimilation_model->get < double >(t, AssimilationModel::INTERC);
        _pari = _assimilation_model->get < double >(t, AssimilationModel::PARI);

        //CulmStockModel
        if(/*_plant_phase != plant::INITIAL and _plant_phase != plant::VEGETATIVE and*/ (_plant_state & plant::INDIV)) {
            _tmp_culm_stock_sum = 0;
            _tmp_culm_deficit_sum = 0;
            _tmp_culm_surplus_sum = 0;
            _tmp_internode_stock_sum = 0;
            _tmp_mainstem_stock_IN = 0;
            _tmp_mainstem_stock = 0;
            _culm_stock_sum = 0;
            _culm_deficit_sum = 0;
            _culm_surplus_sum = 0;
            _internode_stock_sum = 0;
            _mainstem_stock_IN = 0;
            _mainstem_stock = 0;
            _plant_supply = _assimilation_model->get < double >(t, AssimilationModel::ASSIM) + _realloc_sum_supply;
            it = _culm_models.begin();
            while(it != _culm_models.end()) {
                (*it)->stock_model()->put < double >(t, CulmStockModelNG::PLANT_SURPLUS, _stock_model->get < double >(t-1, PlantStockModel::SURPLUS));
                (*it)->stock_model()->put < double >(t, CulmStockModelNG::PLANT_SUPPLY, _plant_supply);
                (*it)->stock_model()->put < double >(t, CulmStockModelNG::PLANT_LEAF_BIOMASS, _leaf_biomass_sum);
                (*it)->compute_stock(t);
                _plant_supply = (*it)->stock_model()->get < double >(t, CulmStockModelNG::NEW_PLANT_SUPPLY);
                _tmp_culm_stock_sum += (*it)->stock_model()->get < double >(t, CulmStockModelNG::CULM_STOCK);
                _tmp_culm_deficit_sum += (*it)->stock_model()->get < double >(t, CulmStockModelNG::CULM_DEFICIT);
                _tmp_internode_stock_sum += (*it)->stock_model()->get< double >(t, CulmStockModelNG::INTERNODE_STOCK);
                if(it == _culm_models.begin()) {
                    _tmp_mainstem_stock = (*it)->stock_model()->get< double >(t, CulmStockModelNG::CULM_STOCK);
                    _tmp_mainstem_stock_IN = (*it)->stock_model()->get< double >(t, CulmStockModelNG::INTERNODE_STOCK);
                }
                ++it;
            }

            it = _culm_models.begin();
            if(_plant_supply > 0) {
                while(it != _culm_models.end()) {
                    if(!((*it)->get < bool, CulmModel >(t, CulmModel::KILL_CULM))) {
                        (*it)->stock_model()->put < double >(t, CulmStockModelNG::PLANT_SUPPLY, _plant_supply);
                        (*it)->stock_model()->iterate_stock(t);
                        _plant_supply = (*it)->stock_model()->get < double >(t, CulmStockModelNG::NEW_PLANT_SUPPLY);
                        _culm_stock_sum += (*it)->stock_model()->get < double >(t, CulmStockModelNG::CULM_STOCK);
                        _culm_deficit_sum += (*it)->stock_model()->get < double >(t, CulmStockModelNG::CULM_DEFICIT);
                        _internode_stock_sum += (*it)->stock_model()->get< double >(t, CulmStockModelNG::INTERNODE_STOCK);
                        if(it == _culm_models.begin()) {
                            _mainstem_stock = (*it)->stock_model()->get< double >(t, CulmStockModelNG::CULM_STOCK);
                            _mainstem_stock_IN = (*it)->stock_model()->get< double >(t, CulmStockModelNG::INTERNODE_STOCK);
                        }
                    }
                    ++it;
                }
            } else {
                _culm_stock_sum = _tmp_culm_stock_sum;
                _culm_deficit_sum = _tmp_culm_deficit_sum;
                _internode_stock_sum = _tmp_internode_stock_sum;
                _mainstem_stock_IN = _tmp_mainstem_stock_IN;
                _mainstem_stock = _tmp_mainstem_stock;
            }
            _culm_surplus_sum = _plant_supply;
        }

        //Root
        _root_model->put < double >(t, RootModel::LEAF_DEMAND_SUM, _leaf_demand_sum);
        _root_model->put < double >(t, RootModel::LEAF_LAST_DEMAND_SUM, _leaf_last_demand_sum);
        _root_model->put < double >(t, RootModel::INTERNODE_DEMAND_SUM, _internode_demand_sum);
        _root_model->put < double >(t, RootModel::INTERNODE_LAST_DEMAND_SUM, _internode_last_demand_sum);
        _root_model->put < plant::plant_state >(t, RootModel::PLANT_STATE, _plant_state);
        _root_model->put < plant::plant_phase >(t, RootModel::PLANT_PHASE, _plant_phase);
        _root_model->put < double >(t, RootModel::CULM_SURPLUS_SUM, _culm_surplus_sum);
        (*_root_model)(t);

        // Stock
        double demand_sum;
        if(_plant_phase == plant::VEGETATIVE) {
            demand_sum = _leaf_demand_sum + _internode_demand_sum + _panicle_demand_sum + _peduncle_demand_sum + _root_model->get < double >(t, RootModel::ROOT_DEMAND);
        } else if(_plant_phase == plant::ELONG or _plant_phase == plant::PI) {
            demand_sum = _leaf_demand_sum + _internode_demand_sum + _panicle_demand_sum + _peduncle_demand_sum;// + _root_model->get < double >(t, RootModel::LAST_ROOT_DEMAND);
        } else {
            demand_sum = _leaf_demand_sum + _internode_demand_sum + _panicle_demand_sum + _peduncle_demand_sum;
        }
        _stock_model->put < double >(t, PlantStockModel::DEMAND_SUM, demand_sum);
        _stock_model->put < double >(t, PlantStockModel::BOOL_CROSSED_PLASTO, _bool_crossed_plasto);
        _stock_model->put < double >(t, PlantStockModel::LEAF_LAST_DEMAND_SUM, _leaf_last_demand_sum);
        _stock_model->put < double >(t, PlantStockModel::PEDUNCLE_LAST_DEMAND_SUM, _peduncle_last_demand_sum);
        _stock_model->put < double >(t, PlantStockModel::INTERNODE_LAST_DEMAND_SUM, _internode_last_demand_sum);
        _stock_model->put < plant::plant_state >(t, PlantStockModel::PLANT_STATE, _plant_state);
        _stock_model->put < double >(t, PlantStockModel::LEAF_BIOMASS_SUM, _leaf_biomass_sum);
        _stock_model->put < double >(t, PlantStockModel::DELETED_LEAF_BIOMASS, _deleted_leaf_biomass);
        _stock_model->put < double >(t, PlantStockModel::REALLOC_BIOMASS_SUM, _realloc_biomass_sum);
        _stock_model->put < double >(t, PlantStockModel::ASSIM,
                                     _assimilation_model->get < double >(t, AssimilationModel::ASSIM));
        _stock_model->put < double >(t, PlantStockModel::CULM_STOCK, _culm_stock_sum);
        _stock_model->put < double >(t, PlantStockModel::CULM_DEFICIT, _culm_deficit_sum);
        if (_plant_phase == plant::PI or _plant_phase == plant::ELONG) {
            _stock_model->put < double >(t, PlantStockModel::CULM_SURPLUS_SUM, _root_model->get< double > (t, RootModel::SURPLUS));
        } else if(_plant_phase == plant::PRE_FLO or _plant_phase == plant::FLO or _plant_phase == plant::END_FILLING or _plant_phase == plant::MATURITY) {
            _stock_model->put < double >(t, PlantStockModel::CULM_SURPLUS_SUM, _culm_surplus_sum);
        } else {
            _stock_model->put < double >(t, PlantStockModel::CULM_SURPLUS_SUM, 0);
        }
        _stock_model->put < plant::plant_phase >(t, PlantStockModel::PLANT_PHASE, _plant_phase);
        (*_stock_model)(t);

        //Variables de visu
        _biomLeaf = _leaf_biomass_sum + _stock_model->get< double > (t, PlantStockModel::STOCK) - _internode_stock_sum;
        _biomin = _internode_biomass_sum + _internode_stock_sum;
        _biominstruct = _internode_biomass_sum;

        // Search leaf to kill
        search_deleted_leaf(t);

        // PHT
        compute_height(t);

        // NB Alive Culms
        double nbc = 0;
        double nbtc = 0;
        _biomAero2 = 0;
        _biomAeroFW = 0;
        _deadleafNb = 0;
        _panicleDW = 0;
        _paniclenb = 0;
        _tillerleafFW = 0;
        std::deque < CulmModel* >::const_iterator itnbc = _culm_models.begin();
        while(itnbc != _culm_models.end()) {
            if(!((*itnbc)->get < bool, CulmModel >(t, CulmModel::KILL_CULM)) and (*itnbc)->get < bool, CulmModel >(t, CulmModel::IS_COMPUTED)) {
                nbc++;
            }
            nbtc++;
            _biomAero2 += (*itnbc)->get < double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM) +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::INTERNODE_BIOMASS_SUM) +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_BIOMASS) +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::PANICLE_WEIGHT);
            _biomAeroFW += (*itnbc)->get < double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM) * _leaf_FW_DW +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::INTERNODE_BIOMASS_SUM) * _internode_FW_DW +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_BIOMASS) * _internode_FW_DW +
                    (*itnbc)->get < double, CulmModel >(t, CulmModel::PANICLE_WEIGHT); //* _panicle_FW_DW;
            _tillerleafFW += (*itnbc)->get < double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM) * _leaf_FW_DW;
            if((*itnbc)->get < bool, CulmModel >(t, CulmModel::CULM_SURVIVED)) {
                _deadleafNb += (*itnbc)->get_dead_phytomer_number(t);
            }
            _panicleDW += (*itnbc)->get < double, CulmModel >(t, CulmModel::PANICLE_WEIGHT);
            _paniclenb += (*itnbc)->get < double, CulmModel >(t, CulmModel::GRAIN_NB) > 0 ? 1 : 0;
            itnbc++;
        }
        _biomAero2 = _biomAero2 +  _stock_model->get< double >(t, PlantStockModel::STOCK);
        _biomAeroFW = _biomAeroFW + (_internode_stock_sum * _internode_FW_DW) + ((_stock_model->get< double >(t, PlantStockModel::STOCK) - _internode_stock_sum) * _leaf_FW_DW);
        _tillerleafFW = _tillerleafFW + ((_stock_model->get< double >(t, PlantStockModel::STOCK) - _internode_stock_sum) * _leaf_FW_DW) - (_biomLeafMainstem * _leaf_FW_DW);
        _biomAeroTot = _biomAero2 + _senesc_dw_sum;

        // VISU
        _tillerNb_1 = nbc;
        _createdTillers = nbtc;
        std::deque < CulmModel* >::const_iterator visumainstem = _culm_models.begin();
        _ms_leaf2_len = (*visumainstem)->get< double, CulmModel >(t, CulmModel::FIRST_LEAF_TOT_LEN); //feuille numéro 2 pour avoir la croissance totale
        _ms_index = (*visumainstem)->get_phytomer_number();
        _biomLeafMainstemstruct = (*visumainstem)->get< double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM);
        _biomLeafMainstem = (*visumainstem)->get< double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM) + (_mainstem_stock - _mainstem_stock_IN);
        _biomInMainstemstruct = (*visumainstem)->get< double, CulmModel >(t, CulmModel::INTERNODE_BIOMASS_SUM);
        _biomInMainstem = (*visumainstem)->get< double, CulmModel >(t, CulmModel::INTERNODE_BIOMASS_SUM) + _mainstem_stock_IN;
        _areaLFEL = (*visumainstem)->get < double, CulmModel >(t, CulmModel::LAST_LEAF_BLADE_AREA);
        _nbleaf = (*visumainstem)->get < double, CulmModel >(t, CulmModel::NBLEAF); // nombre de feuilles vertes (>= 50% blade area)
        _internode_length_mainstem = (*visumainstem)->get < double, CulmModel >(t, CulmModel::INTERNODE_LEN_SUM);
        _total_length_mainstem = (*visumainstem)->get < double, CulmModel >(t, CulmModel::INTERNODE_LEN_SUM) + (*visumainstem)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_LEN);
        _panicleMainstemDW = (*visumainstem)->get < double, CulmModel >(t, CulmModel::PANICLE_WEIGHT);
        _biomMainstem = _biomLeafMainstem + _biomInMainstem + _panicleMainstemDW;
        _biomMainstemFW = (_biomLeafMainstem * _leaf_FW_DW) + (_biomInMainstem * _internode_FW_DW); //+ _panicleMainstemDW * _panicle_FW_DW;
        _biomBladeMainstemFW = _biomLeafMainstem * _G_L * _leaf_FW_DW;
        _tillerFW = _biomAeroFW - _biomMainstemFW;
        if(_biomLeaf > 0) {
            _slaplant = _leaf_blade_area_sum / (_biomLeaf * _G_L);
        } else {
            _slaplant = 0;
        }

        _biomLeafTot = _biomLeaf + _senesc_dw_sum;
        _biomInSheathMainstem = _biomLeafMainstem - (_biomLeafMainstem * _G_L) + _biomInMainstem;
        _biomInSheath = _biomLeaf - (_biomLeaf * _G_L) + _biomin;
    }

    void create_culm(double t, int n)
    {
        for (int i = 0; i < n; ++i) {
            CulmModel* meristem = new CulmModel(_culm_models.size() + 1);
            setsubmodel(CULMS, meristem);
            meristem->put(t, CulmModel::LL_BL, _LL_BL);
            meristem->put(t, CulmModel::PLANT_PHENOSTAGE, _phenostage);
            meristem->init(t, _parameters);
            _culm_models.push_back(meristem);
        }
    }

    void compute_culms(double t)
    {
        std::deque < CulmModel* >::const_iterator it = _culm_models.begin();
        while (it != _culm_models.end()) {
            (*it)->put(t, CulmModel::MS_PHYT_INDEX, _ms_index);
            (*it)->put(t, CulmModel::PLANT_BOOL_CROSSED_PLASTO, _bool_crossed_plasto);
            (*it)->put(t, CulmModel::MS_SHEATH_LLL, _sheath_LLL);
            (*it)->put(t, CulmModel::DD, _DD);
            (*it)->put(t, CulmModel::EDD, _EDD);
            (*it)->put(t, CulmModel::DELTA_T, _deltaT);
            (*it)->put(t, CulmModel::FTSW, _water_balance_model->get < double >(t, WaterBalanceModel::FTSW));
            (*it)->put(t, CulmModel::FCSTR, _water_balance_model->get < double >(t, WaterBalanceModel::FCSTR));
            (*it)->put(t, CulmModel::FCSTRI, _water_balance_model->get < double >(t, WaterBalanceModel::FCSTRI));
            (*it)->put(t, CulmModel::FCSTRL, _water_balance_model->get < double >(t, WaterBalanceModel::FCSTRL));
            (*it)->put(t, CulmModel::FCSTRLLEN, _water_balance_model->get < double >(t, WaterBalanceModel::FCSTRLLEN));
            (*it)->put < int > (t, CulmModel::PLANT_PHENOSTAGE, _phenostage);
            (*it)->put < int > (t, CulmModel::PLANT_APPSTAGE, _appstage);
            (*it)->put < int > (t, CulmModel::PLANT_LIGSTAGE, _ligstage);
            (*it)->put(t, CulmModel::PREDIM_LEAF_ON_MAINSTEM, _predim_leaf_on_mainstem);
            (*it)->put(t, CulmModel::PREDIM_APP_LEAF_MS, _predim_app_leaf_on_mainstem);
            (*it)->put(t, CulmModel::SLA, _sla);
            (*it)->put < plant::plant_state >(t, CulmModel::PLANT_STATE, _plant_state);
            (*it)->put < plant::plant_phase >(t, CulmModel::PLANT_PHASE, _plant_phase);
            (*it)->put(t, CulmModel::TEST_IC, _stock_model->get < double >(t-1, PlantStockModel::TEST_IC));
            (*it)->put(t, CulmModel::PLANT_STOCK, _stock);
            (*it)->put(t, CulmModel::PLANT_DEFICIT, _deficit);
            (*it)->put(t, CulmModel::ASSIM, _assimilation_model->get < double >(t-1, AssimilationModel::ASSIM));
            (*it)->put(t, CulmModel::MGR, _MGR);
            (*it)->put(t, CulmModel::LL_BL, _LL_BL);
            (*it)->put(t, CulmModel::IS_FIRST_DAY_PI, _is_first_day_pi);
            (**it)(t);
            ++it;
        }

        _leaf_biomass_sum = 0;
        _last_leaf_biomass_sum = 0;
        _leaf_last_demand_sum = 0;
        _leaf_demand_sum = 0;
        _internode_last_demand_sum = 0;
        _internode_demand_sum = 0;
        _internode_biomass_sum = 0;
        _leaf_blade_area_sum = 0;
        _realloc_biomass_sum = 0;
        _senesc_dw_sum = 0;
        _senesc_dw = 0;
        _panicle_demand_sum = 0;
        _peduncle_demand_sum = 0;
        _peduncle_last_demand_sum = 0;
        _peduncle_biomass_sum = 0;
        _realloc_sum_supply = 0;

        it = _culm_models.begin();
        _predim_leaf_on_mainstem = (*it)->get < double, CulmModel > (t, CulmModel::STEM_LEAF_PREDIM);
        _predim_app_leaf_on_mainstem = (*it)->get < double, CulmModel > (t, CulmModel::STEM_APP_LEAF_PREDIM);
        _sheath_LLL = (*it)->get < double, CulmModel >(t, CulmModel::SHEATH_LLL);
        _lig_index = (*it)->get < double, CulmModel >(t, CulmModel::LIG_INDEX);

        while (it != _culm_models.end()) {
            _panicle_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::PANICLE_DAY_DEMAND);
            _peduncle_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_DAY_DEMAND);
            _peduncle_last_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_LAST_DEMAND);
            _leaf_biomass_sum += (*it)->get < double, CulmModel >(t, CulmModel::LEAF_BIOMASS_SUM);
            _last_leaf_biomass_sum += (*it)->get < double, CulmModel >(t, CulmModel::LAST_LEAF_BIOMASS_SUM);
            _leaf_last_demand_sum += (*it)->get < double, CulmModel>(t, CulmModel::LEAF_LAST_DEMAND_SUM);
            _leaf_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::LEAF_DEMAND_SUM);
            _internode_last_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::INTERNODE_LAST_DEMAND_SUM);
            _internode_demand_sum += (*it)->get < double, CulmModel >(t, CulmModel::INTERNODE_DEMAND_SUM);
            _internode_biomass_sum += (*it)->get < double, CulmModel >(t, CulmModel::INTERNODE_BIOMASS_SUM);
            _peduncle_biomass_sum += (*it)->get < double, CulmModel >(t, CulmModel::PEDUNCLE_BIOMASS);
            _leaf_blade_area_sum += (*it)->get < double, CulmModel>(t, CulmModel::LEAF_BLADE_AREA_SUM);
            _realloc_biomass_sum += (*it)->get < double, CulmModel>(t, CulmModel::REALLOC_BIOMASS_SUM);
            _senesc_dw_sum += (*it)->get < double, CulmModel>(t, CulmModel::SENESC_DW_SUM) + (*it)->get < double, CulmModel>(t, CulmModel::DELETED_SENESC_DW_SUM);
            _senesc_dw += (*it)->get < double, CulmModel>(t, CulmModel::SENESC_DW) + (*it)->get < double, CulmModel>(t, CulmModel::DELETED_SENESC_DW);
            _realloc_sum_supply += (*it)->get < double, CulmModel >(t, CulmModel::REALLOC_SUPPLY);
            _leaf_delay = (*it)->get< double, CulmModel>(t, CulmModel::LEAF_DELAY);
            ++it;
        }
    }

    void compute_height(double t)
    {
        _height = 0;
        if (_culm_models.empty()) {
            return;
        }
        auto it = _culm_models.begin();
        _height += (*it)->get < double, CulmModel >(t, CulmModel::INTERNODE_LEN_SUM);
        _height_ped = _height;
        if ((*it)->get_phytomer_number() == 1) {
            _height += (*it)->get < double, CulmModel >(t, CulmModel::FIRST_LEAF_LEN);
            _height_ped = _height;
        } else {
            double tmp = (*it)->get < double, CulmModel >(t, CulmModel::LAST_LIGULATED_LEAF_SHEATH_LEN);
            _height += (*it)->get < double, CulmModel >(t, CulmModel::LAST_LIGULATED_LEAF_SHEATH_LEN);
            if (tmp > (*it)->get< double, CulmModel >(t, CulmModel::PEDUNCLE_LEN)) {
                _height_ped += tmp;
            } else {
                _height_ped += (*it)->get< double, CulmModel >(t, CulmModel::PEDUNCLE_LEN);
            }
        }
    }

    void delete_leaf(double t)
    {
        if (_culm_index != -1 and _leaf_index != -1) {
            _culm_models[_culm_index]->delete_leaf(t, _leaf_index, _deleted_leaf_biomass, _deleted_internode_biomass);
            _leaf_blade_area_sum -= _deleted_leaf_blade_area;
        }
    }

    void search_deleted_leaf(double t)
    {
        _culm_index = -1;
        _leaf_index = -1;
        _deleted_leaf_biomass = 0;
        _deleted_internode_biomass = 0;
        _deleted_leaf_blade_area = 0;
        _qty = 0;
        if (_stock_model->get < double >(t, PlantStockModel::STOCK) == 0) {
            std::deque < CulmModel* >::const_iterator it = _culm_models.begin();
            if(/*_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE*/ !(_plant_state & plant::INDIV)) {
                double tmp_date = t;
                std::deque < CulmModel* >::const_iterator it = _culm_models.begin();
                int i = 0;
                while (it != _culm_models.end()) {
                    double creation_date = (*it)->get_first_alive_leaf_creation_date(t);
                    if(creation_date < tmp_date) {
                        _culm_index = i;
                        _leaf_index = (*it)->get_first_alive_leaf_index(t);
                        tmp_date = creation_date;
                    }
                    ++it;
                    ++i;
                }
            }
            if (_leaf_index != -1) {
                _deleted_leaf_biomass =
                        _culm_models[_culm_index]->get_leaf_biomass(t, _leaf_index);
                _deleted_leaf_blade_area =
                        _culm_models[_culm_index]->get_leaf_blade_area(t,_leaf_index);
                if (_culm_models[_culm_index]->get_alive_phytomer_number() < 2 and _culm_index == 0) {
                    _plant_state << plant::NOGROWTH;
                    _stock = 0;
                }
            }
        }
    }

    void init(double t, const ecomeristem::ModelParameters& parameters)
    {
        //parameters
        _parameters = parameters;
        _nbleaf_enabling_tillering = _parameters.get("nb_leaf_enabling_tillering");
        _LL_BL_init = _parameters.get("LL_BL_init");
        _nb_leaf_param2 = _parameters.get("nb_leaf_param2");
        _slope_LL_BL_at_PI = _parameters.get("slope_LL_BL_at_PI");
        _nb_leaf_enabling_tillering = _parameters.get("nb_leaf_enabling_tillering");
        _coeff_MGR_PI = _parameters.get("coef_MGR_PI");
        _nb_leaf_stem_elong = _parameters.get("nb_leaf_stem_elong");
        _phenostage_pre_flo_to_flo  = _parameters.get("phenostage_PRE_FLO_to_FLO");
        _phenostage_to_end_filling = _parameters.get("phenostage_to_end_filling");
        _phenostage_to_maturity = _parameters.get("phenostage_to_maturity");
        _Ict = _parameters.get("Ict");
        _leaf_stock_max = _parameters.get("leaf_stock_max");
        _realocationCoeff = _parameters.get("realocationCoeff");
        _FSLA = _parameters.get("FSLA");
        _plasto_init = _parameters.get("plasto_init");
        _ligulo_init = _parameters.get("ligulo_init");
        _phyllo_init = _parameters.get("phyllo_init");
        _SLAp = _parameters.get("SLAp");
        _Tb = _parameters.get("Tb");
        _maxleaves = _parameters.get("maxleaves");
        _G_L = parameters.get("G_L");
        _intercmodel = parameters.get("intercmodel");
        _leaf_FW_DW = parameters.get("leaf_FW_DW");
        _internode_FW_DW = parameters.get("internode_FW_DW");
        _nb_leaf_indiv = parameters.get("nb_leaf_indiv");

        //Attributes for culmmodel
        _LL_BL = _LL_BL_init;
        _phenostage = 4;

        //local init
        CulmModel* meristem = new CulmModel(1);
        setsubmodel(CULMS, meristem);
        meristem->put(t, CulmModel::LL_BL, _LL_BL);
        meristem->put(t, CulmModel::PLANT_PHENOSTAGE, _phenostage);
        meristem->init(t, parameters);
        _culm_models.push_back(meristem);

        //submodels
        _water_balance_model->init(t, parameters);
        _stock_model->init(t, parameters);
        _assimilation_model->init(t, parameters);
        _interception_model->init(t, parameters);
        _root_model->init(t, parameters);

        //vars
        _predim_leaf_on_mainstem = 0;

        //internal variables (local)
        _tae = 0;
        _nb_tillers = 0;
        _paniclenb = 0;
        _nbExistingTillers = 1;
        _lig = 0;
        _app = 1;
        _lig_1 = 0;
        _leaf_biomass_sum = 0;
        _leaf_demand_sum = 0;
        _leaf_last_demand_sum = 0;
        _internode_demand_sum = 0;
        _internode_last_demand_sum = 0;
        _internode_biomass_sum = 0;
        _leaf_blade_area_sum = 0;
        _senesc_dw_sum = 0;
        _culm_stock_sum = 0;
        _culm_deficit_sum = 0;
        _culm_surplus_sum = 0;
        _panicle_demand_sum = 0;
        _plant_phase = plant::VEGETATIVE;
        _plant_state = plant::NO_STATE;
        _height = 0;
        _height_ped = 0;
        _MGR = parameters.get("MGR_init");
        _TT_lig = 0;
        _IH = 0;
        _biomLeaf = 0;
        _last_leaf_biomass_sum = 0;
        _is_first_day_pi = false;
        _internode_stock_sum = 0;
        _peduncle_biomass_sum = 0;
        _peduncle_demand_sum = 0;
        _peduncle_last_demand_sum = 0;
        _biomin = 0;
        _biominstruct = 0;
        _deleted_leaf_biomass = 0;
        _deleted_leaf_blade_area = 0;
        _culm_index = -1;
        _leaf_index = -1;
        _stock = 0;
        _deficit = 0;
        _qty = 0;
        _deleted_internode_biomass = 0;
        _last_time = 0;
        _tmp_culm_stock_sum = 0;
        _tmp_culm_deficit_sum = 0;
        _tmp_culm_surplus_sum = 0;
        _tmp_internode_stock_sum = 0;
        _plant_supply = 0;
        _realloc_sum_supply = 0;
        _Ta = 0;
        _deltaT = 0;
        _TT = 0;
        _bool_crossed_plasto = 0;
        _bool_crossed_phyllo = 0;
        _bool_crossed_ligulo = 0;
        _EDD = 0;
        _DD = 0;
        _ligulo_visu = _ligulo_init;
        _plasto_visu = _plasto_init;
        _phyllo_visu = _phyllo_init;
        _phenostage = 4;
        _phenostage_cste = _phenostage;
        _sla = _FSLA;
        _tillerNb_1 = 1;
        _nbleaf = 1;
        _biomAero2 = 0;
        _biomLeafMainstem = 0;
        _biomInMainstem = 0;
        _areaLFEL = 0;
        _mainstem_stock_IN = 0;
        _tmp_mainstem_stock_IN = 0;
        _biomInMainstemstruct = 0;
        _biomLeafMainstemstruct = 0;
        _mainstem_stock = 0;
        _tmp_mainstem_stock = 0;
        _deadleafNb = 0;
        _internode_length_mainstem = 0;
        _panicleMainstemDW = 0;
        _panicleDW = 0;
        _leaf_delay = 0;
        _DD_phyllo = 0;
        _DD_ligulo = 0;
        _sheath_LLL = 0;
        _phenostage_at_flo = 0;
        _lig_index = 0;
        _realloc_biomass_sum = 0;
        _appstage = 1;
        _ligstage = 0;
        _visi = 0;
        _ms_index = 1;
        _predim_app_leaf_on_mainstem = 0;
        _total_length_mainstem = 0;
        _biomMainstem = 0;
        _slaplant = 0;
        _biomLeafTot = 0;
        _senesc_dw = 0;
        _biomInSheathMainstem = 0;
        _biomInSheath = 0;
        _biomAeroTot = 0;
        _ms_leaf2_len = 0;
        _interc1 = 0;
        _interc2 = 0;
        _pari = 0;
        _biomAeroFW = 0;
        _biomMainstemFW = 0;
        _tillerFW = 0;
        _biomBladeMainstemFW = 0;
        _tillerleafFW = 0;
        _is_first_day_passed = false;
        _createdTillers = 1;
        //_current_date = t - parameters.beginDate;
    }

private:
    double _last_time;

    ecomeristem::ModelParameters _parameters;
    // submodels
    std::deque < CulmModel* > _culm_models;
    std::unique_ptr < model::WaterBalanceModel > _water_balance_model;
    std::unique_ptr < model::PlantStockModel > _stock_model;
    std::unique_ptr < model::AssimilationModel > _assimilation_model;
    std::unique_ptr < model::InterceptionModel > _interception_model;
    std::unique_ptr < model::RootModel > _root_model;

    // parameters
    double _coeff_MGR_PI;
    double _nbleaf_enabling_tillering;
    double _nb_leaf_param2;
    double _slope_LL_BL_at_PI;
    double _LL_BL_init;
    double _nb_leaf_stem_elong;
    double _phenostage_pre_flo_to_flo;
    double _phenostage_to_end_filling;
    double _phenostage_to_maturity;
    double _Ict;
    double _leaf_stock_max;
    double _nb_leaf_enabling_tillering;
    double _realocationCoeff;
    double _FSLA;
    double _plasto_init;
    double _ligulo_init;
    double _phyllo_init;
    double _SLAp;
    double _Tb;
    double _maxleaves;
    double _tae;
    double _G_L;
    double _intercmodel;
    double _leaf_FW_DW;
    double _internode_FW_DW;
    double _nb_leaf_indiv;


    // vars
    double _predim_leaf_on_mainstem;
    // internals
    double _nb_tillers;
    double _nbExistingTillers;
    double _MGR;
    double _LL_BL;
    double _lig;
    double _app;
    double _leaf_biomass_sum;
    double _leaf_demand_sum;
    double _leaf_blade_area_sum;
    double _leaf_last_demand_sum;
    double _internode_biomass_sum;
    double _internode_demand_sum;
    double _internode_last_demand_sum;
    double _senesc_dw_sum;
    double _realloc_biomass_sum;
    double _culm_stock_sum;
    double _culm_deficit_sum;
    double _culm_surplus_sum;
    double _panicle_demand_sum;
    double _height;
    double _height_ped;
    double _lig_1;
    double _TT_lig;
    double _IH;
    double _biomLeaf;
    double _last_leaf_biomass_sum;
    bool _is_first_day_pi;
    double _internode_stock_sum;
    double _biomin;
    double _biominstruct;
    double _peduncle_biomass_sum;
    double _peduncle_demand_sum;
    double _peduncle_last_demand_sum;
    double _deleted_leaf_biomass;
    double _deleted_leaf_blade_area;
    int _culm_index;
    int _leaf_index;
    double _stock;
    double _deficit;
    double _qty;
    double _deleted_internode_biomass;
    double _tmp_culm_stock_sum;
    double _tmp_culm_deficit_sum;
    double _tmp_culm_surplus_sum;
    double _tmp_internode_stock_sum;
    double _plant_supply;
    double _realloc_sum_supply;
    double _tillerNb_1;
    double _nbleaf;
    double _biomAero2;
    double _biomLeafMainstem;
    double _biomInMainstem;
    double _areaLFEL;
    double _mainstem_stock_IN;
    double _tmp_mainstem_stock_IN;
    double _biomInMainstemstruct;
    double _Ta;
    double _deltaT;
    double _TT;
    double _bool_crossed_plasto;
    double _bool_crossed_phyllo;
    double _bool_crossed_ligulo;
    double _EDD;
    double _DD;
    double _ligulo_visu;
    double _plasto_visu;
    int _phenostage;
    double _sla;
    double _phenostage_cste;
    double _mainstem_stock;
    double _tmp_mainstem_stock;
    double _biomLeafMainstemstruct;
    int _deadleafNb;
    double _internode_length_mainstem;
    double _panicleMainstemDW;
    double _panicleDW;
    double _leaf_delay;
    double _DD_phyllo;
    double _DD_ligulo;
    double _phyllo_visu;
    double _sheath_LLL;
    double _appstage;
    double _appstage_cste;
    double _ligstage;
    double _ligstage_cste;
    double _phenostage_at_flo;
    double _lig_index;
    int _ms_index;
    double _visi;
    double _predim_app_leaf_on_mainstem;
    double _paniclenb;
    double _total_length_mainstem;
    double _biomMainstem;
    double _slaplant;
    double _biomLeafTot;
    double _senesc_dw;
    double _biomInSheathMainstem;
    double _biomInSheath;
    double _biomAeroTot;
    double _ms_leaf2_len;
    double _pari;
    double _biomAeroFW;
    double _biomMainstemFW;
    double _tillerFW;
    double _biomBladeMainstemFW;
    double _tillerleafFW;
    bool _is_first_day_passed;
    double _createdTillers;
    //double _current_date;

    // test
    double _interc1;
    double _interc2;

    //internal states
    plant::plant_state _plant_state;
    plant::plant_phase _plant_phase;
};

#endif //PLANT_MODEL_HPP
