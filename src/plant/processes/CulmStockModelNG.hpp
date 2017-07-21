#ifndef CULMSTOCKMODELNG_H
#define CULMSTOCKMODELNG_H

#include <defines.hpp>

namespace model {

class CulmStockModelNG : public AtomicModel < CulmStockModelNG >
{
public:
    enum internals { MAX_RESERVOIR_DISPO_INTERNODE, RESERVOIR_DISPO_INTERNODE,
                     MAX_RESERVOIR_DISPO_LEAF, RESERVOIR_DISPO_LEAF,
                     INTERNODE_STOCK, LEAF_STOCK, DEMAND_INTERNODE_STORAGE,
                     REMAIN_TO_STORE, REMAIN_TO_STORE_1, CULM_DEMAND_SUM,
                     CULM_SUPPLY, CULM_IC, CULM_DEFICIT, CULM_DEFICIT_1,
                     INTERMEDIATE, INTERMEDIATE2, INTERMEDIATE3,
                     NEW_PLANT_SUPPLY, CULM_STOCK, CULM_DEMAND, LEAF_STOCK_INIT,
                     CULM_SURPLUS, MAX_RESERVOIR_DISPO };


    enum externals { LEAF_DEMAND_SUM, INTERNODE_DEMAND_SUM,
                     PANICLE_DEMAND, PEDUNCLE_DEMAND, LAST_DEMAND, REALLOC_BIOMASS,
                     PLANT_LEAF_BIOMASS, PLANT_STOCK, PLANT_SURPLUS, PLANT_SUPPLY,
                     INTERNODE_BIOMASS_SUM, LEAF_BIOMASS_SUM, PLANT_PHASE,
                     IS_FIRST_DAY_OF_INDIVIDUALIZATION, KILL_CULM };


    CulmStockModelNG() {
        Internal(MAX_RESERVOIR_DISPO_INTERNODE, &CulmStockModelNG::_max_reservoir_dispo_internode);
        Internal(RESERVOIR_DISPO_INTERNODE, &CulmStockModelNG::_reservoir_dispo_internode);
        Internal(MAX_RESERVOIR_DISPO_LEAF, &CulmStockModelNG::_max_reservoir_dispo_leaf);
        Internal(RESERVOIR_DISPO_LEAF, &CulmStockModelNG::_reservoir_dispo_leaf);
        Internal(INTERNODE_STOCK, &CulmStockModelNG::_internode_stock);
        Internal(LEAF_STOCK, &CulmStockModelNG::_leaf_stock);
        Internal(DEMAND_INTERNODE_STORAGE, &CulmStockModelNG::_demand_internode_storage);
        Internal(REMAIN_TO_STORE, &CulmStockModelNG::_remain_to_store);
        Internal(REMAIN_TO_STORE_1, &CulmStockModelNG::_remain_to_store_1);
        Internal(CULM_DEMAND_SUM, &CulmStockModelNG::_culm_demand_sum);
        Internal(CULM_SUPPLY, &CulmStockModelNG::_culm_supply);
        Internal(CULM_IC, &CulmStockModelNG::_culm_ic);
        Internal(CULM_DEFICIT, &CulmStockModelNG::_culm_deficit);
        Internal(CULM_DEFICIT_1, &CulmStockModelNG::_culm_deficit_1);
        Internal(INTERMEDIATE, &CulmStockModelNG::_intermediate);
        Internal(INTERMEDIATE2, &CulmStockModelNG::_intermediate2);
        Internal(INTERMEDIATE3, &CulmStockModelNG::_intermediate3);
        Internal(NEW_PLANT_SUPPLY, &CulmStockModelNG::_new_plant_supply);
        Internal(CULM_STOCK, &CulmStockModelNG::_culm_stock);
        Internal(CULM_DEMAND, &CulmStockModelNG::_culm_demand);
        Internal(LEAF_STOCK_INIT, &CulmStockModelNG::_leaf_stock_init);
        Internal(CULM_SURPLUS, &CulmStockModelNG::_culm_surplus);
        Internal(MAX_RESERVOIR_DISPO, &CulmStockModelNG::_max_reservoir_dispo);

        External(LEAF_DEMAND_SUM, &CulmStockModelNG::_leaf_demand_sum);
        External(INTERNODE_DEMAND_SUM, &CulmStockModelNG::_internode_demand_sum);
        External(PANICLE_DEMAND, &CulmStockModelNG::_panicle_demand);
        External(PEDUNCLE_DEMAND, &CulmStockModelNG::_peduncle_demand);
        External(REALLOC_BIOMASS, &CulmStockModelNG::_realloc_biomass);
        External(PLANT_LEAF_BIOMASS, &CulmStockModelNG::_plant_leaf_biomass);
        External(LAST_DEMAND, &CulmStockModelNG::_last_demand);
        External(PLANT_STOCK, &CulmStockModelNG::_plant_stock);
        External(PLANT_SURPLUS, &CulmStockModelNG::_plant_surplus);
        External(PLANT_SUPPLY, &CulmStockModelNG::_plant_supply);
        External(INTERNODE_BIOMASS_SUM, &CulmStockModelNG::_internode_biomass_sum);
        External(LEAF_BIOMASS_SUM, &CulmStockModelNG::_leaf_biomass_sum);
        External(PLANT_PHASE, &CulmStockModelNG::_plant_phase);
        External(IS_FIRST_DAY_OF_INDIVIDUALIZATION, &CulmStockModelNG::_is_first_day_of_individualization);
        External(KILL_CULM, &CulmStockModelNG::_kill_culm);

    }

    virtual ~CulmStockModelNG()
    {}


    void compute(double t, bool /* update */) {
        if (_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE or _kill_culm) {
            _max_reservoir_dispo_internode = 0;
            _reservoir_dispo_internode = 0;
            _max_reservoir_dispo_leaf = 0;
            _reservoir_dispo_leaf = 0;
            _internode_stock = 0;
            _leaf_stock = 0;
            _demand_internode_storage = 0;
            _remain_to_store = 0;
            _remain_to_store_1 = 0;
            _culm_demand_sum = 0;
            _culm_supply = 0;
            _culm_ic = 0;
            _culm_deficit = 0;
            _culm_deficit_1 = 0;
            _culm_surplus = 0;
            _intermediate = 0;
            _intermediate2 = 0;
            _intermediate3 = 0;
            _culm_stock = 0;
            _leaf_stock_init = 0;
            _new_plant_supply = _plant_supply;
            return;
        }
        if(_is_first_day_of_individualization) {
            //First day, _leaf_stock initialisation
            _leaf_stock_init = _plant_stock * (_leaf_biomass_sum / _plant_leaf_biomass);
            _leaf_stock = _leaf_stock_init;
        }
        _max_reservoir_dispo_internode = _maximum_reserve_in_internode * _internode_biomass_sum;
        _max_reservoir_dispo_leaf = _leaf_stock_max * _leaf_biomass_sum;
        _max_reservoir_dispo = (_maximum_reserve_in_internode *
                                _internode_biomass_sum) + (_leaf_stock_max * _leaf_biomass_sum);

        _reservoir_dispo_internode = _max_reservoir_dispo_internode - _internode_stock;
        _culm_demand = _leaf_demand_sum + _panicle_demand + _peduncle_demand + _internode_demand_sum + _last_demand;
        _demand_internode_storage = _reservoir_dispo_internode * _coeff_active_storage_IN;
        _culm_demand_sum = _culm_demand + _demand_internode_storage;
        _culm_supply = std::min(_plant_supply, _culm_demand_sum);
        _culm_ic = std::min(_culm_supply / std::max(_culm_demand_sum, 0.00001), 5.);
        _intermediate = _culm_supply - _culm_demand_sum + _realloc_biomass + _leaf_stock;
        _intermediate2 = std::max((1 - _coeff_remob) * _internode_stock, std::min(_internode_stock + _reservoir_dispo_internode, _internode_stock + _intermediate + _culm_deficit + _demand_internode_storage));
        _culm_deficit_1 = _culm_deficit;
        _culm_deficit = std::min(0., _culm_deficit + _intermediate + (_coeff_remob * _internode_stock) + _demand_internode_storage);
        _remain_to_store = std::max(0., _intermediate + _demand_internode_storage - _reservoir_dispo_internode + _culm_deficit_1 - _culm_deficit);
        _remain_to_store_1 = _remain_to_store;
        _internode_stock = _intermediate2;
        _leaf_stock = std::min(_max_reservoir_dispo_leaf, _remain_to_store);
        _remain_to_store = _remain_to_store - _leaf_stock;
        _reservoir_dispo_internode = _max_reservoir_dispo_internode - _internode_stock;
        _reservoir_dispo_leaf = _max_reservoir_dispo_leaf - _leaf_stock;
        _culm_stock = _leaf_stock + _internode_stock;
        _new_plant_supply = _plant_supply - _culm_supply + _remain_to_store;
        _culm_surplus = std::max(0., _culm_stock - _culm_demand + _culm_supply - _max_reservoir_dispo + _realloc_biomass);
    }

    void iterate_stock(double t) {
        _culm_deficit_1 = _culm_deficit;
        _culm_deficit = std::min(0., _culm_deficit + _plant_supply);
        _new_plant_supply = std::max(0.,_plant_supply + _culm_deficit_1);
        _internode_stock = _internode_stock + std::min(_reservoir_dispo_internode, _new_plant_supply);
        _intermediate3 = _new_plant_supply - std::min(_reservoir_dispo_internode, _new_plant_supply);
        _leaf_stock = _leaf_stock + std::min(_reservoir_dispo_leaf, _intermediate3);
        _new_plant_supply = _intermediate3 - std::min(_reservoir_dispo_leaf, _intermediate3);
        _reservoir_dispo_internode = _max_reservoir_dispo_internode - _internode_stock;
        _reservoir_dispo_leaf = _max_reservoir_dispo_leaf - _leaf_stock;
        _culm_stock = _leaf_stock + _internode_stock;
        _culm_surplus = std::max(0., _culm_stock - _culm_demand + _culm_supply - _max_reservoir_dispo + _realloc_biomass);
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        //permet le passage du get à t0 en mimant un isComputed au temps t
        last_time = t-1;

        // parameters
        _parameters = parameters;
        _maximum_reserve_in_internode = parameters.get("maximumReserveInInternode");
        _leaf_stock_max = parameters.get("leaf_stock_max");
        _coeff_remob = parameters.get("coeff_remob");
        _coeff_active_storage_IN = parameters.get("coeff_active_storage_IN");

        // internals            
        _max_reservoir_dispo_internode = 0;
        _reservoir_dispo_internode = 0;
        _max_reservoir_dispo_leaf = 0;
        _reservoir_dispo_leaf = 0;
        _internode_stock = 0;
        _leaf_stock = 0;
        _demand_internode_storage = 0;
        _remain_to_store = 0;
        _remain_to_store_1 = 0;
        _culm_demand_sum = 0;
        _culm_supply = 0;
        _culm_ic = 0;
        _culm_deficit = 0;
        _culm_deficit_1 = 0;
        _intermediate = 0;
        _intermediate2 = 0;
        _intermediate3 = 0;
        _new_plant_supply = 0;
        _culm_stock = 0;
        _leaf_stock_init = 0;
        _culm_surplus = 0;
        _max_reservoir_dispo = 0;

    }

private:
    ecomeristem::ModelParameters _parameters;

    //parameters
    double _coeff_remob;
    double _coeff_active_storage_IN;
    double _nb_leaf_stem_elong;
    double _nb_leaf_pi;
    double _maximum_reserve_in_internode;
    double _leaf_stock_max;

    //Internals
    double _max_reservoir_dispo_internode;
    double _reservoir_dispo_internode;
    double _max_reservoir_dispo_leaf;
    double _reservoir_dispo_leaf;
    double _internode_stock;
    double _leaf_stock;
    double _demand_internode_storage;
    double _remain_to_store;
    double _remain_to_store_1;
    double _culm_demand_sum;
    double _culm_supply;
    double _culm_ic;
    double _culm_deficit;
    double _culm_deficit_1;
    double _culm_stock;
    double _intermediate;
    double _intermediate2;
    double _intermediate3;
    double _new_plant_supply;
    double _culm_demand;
    double _leaf_stock_init;
    double _culm_surplus;
    double _max_reservoir_dispo;

    //Externals
    double _leaf_demand_sum;
    double _internode_demand_sum;
    double _panicle_demand;
    double _peduncle_demand;
    double _last_demand;
    double _realloc_biomass;
    double _plant_leaf_biomass;
    double _plant_stock;
    double _plant_surplus;
    double _plant_supply;
    double _internode_biomass_sum;
    double _leaf_biomass_sum;
    bool _is_first_day_of_individualization;
    bool _kill_culm;
    plant::plant_phase _plant_phase;

};

} // namespace model

#endif // CULMSTOCKMODELNG_H
