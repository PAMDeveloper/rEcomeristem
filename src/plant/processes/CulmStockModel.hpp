#ifndef CULMSTOCKMODEL_H
#define CULMSTOCKMODEL_H

#include <defines.hpp>

namespace model {

class CulmStockModel : public AtomicModel < CulmStockModel >
{
public:
    enum internals { STOCK, SUPPLY, MAX_RESERVOIR_DISPO, INTERMEDIATE, DEFICIT, SURPLUS,
                     FIRST_DAY, LEAF_STOCK, STOCK_CULM, STOCK_INTERNODE,
                     CULM_BIOMASS, DEFICIT_CULM };

    enum externals { ASSIM, LEAF_BIOMASS_SUM, PLANT_LEAF_BIOMASS_SUM,
                     INTERNODE_BIOMASS_SUM, PLANT_BIOMASS_SUM, PLANT_STOCK,
                     PLANT_DEFICIT, INTERNODE_DEMAND_SUM, LEAF_DEMAND_SUM,
                     INTERNODE_LAST_DEMAND_SUM, LEAF_LAST_DEMAND_SUM,
                     REALLOC_BIOMASS_SUM, PLANT_PHASE, CULM_PHASE,
                     PANICLE_DAY_DEMAND, PANICLE_WEIGHT, LAST_LEAF_BIOMASS_SUM,
                     LAST_PLANT_LEAF_BIOMASS_SUM, IS_FIRST_DAY_PI,
                     PEDUNCLE_LAST_DEMAND, PEDUNCLE_DAY_DEMAND};


    CulmStockModel() {
        Internal(STOCK, &CulmStockModel::_stock);
        Internal(SUPPLY, &CulmStockModel::_supply);
        Internal(MAX_RESERVOIR_DISPO, &CulmStockModel::_max_reservoir_dispo);
        Internal(DEFICIT, &CulmStockModel::_deficit);
        Internal(SURPLUS, &CulmStockModel::_surplus);
        Internal(FIRST_DAY, &CulmStockModel::_first_day);
        Internal(INTERMEDIATE, &CulmStockModel::_intermediate);
        Internal(LEAF_STOCK, &CulmStockModel::_leaf_stock);
        Internal(STOCK_CULM, &CulmStockModel::_stock_culm);
        Internal(CULM_BIOMASS, &CulmStockModel::_culm_biomass);
        Internal(STOCK_INTERNODE, &CulmStockModel::_stock_internode);
        Internal(DEFICIT_CULM, &CulmStockModel::_deficit_culm);

        External(PLANT_BIOMASS_SUM, &CulmStockModel::_plant_biomass_sum);
        External(LAST_PLANT_LEAF_BIOMASS_SUM, &CulmStockModel::_last_plant_biomass_sum);
        External(PLANT_LEAF_BIOMASS_SUM, &CulmStockModel::_plant_leaf_biomass_sum);
        External(LAST_LEAF_BIOMASS_SUM, &CulmStockModel::_last_leaf_biomass_sum);
        External(ASSIM, &CulmStockModel::_assim);
        External(LEAF_BIOMASS_SUM, &CulmStockModel::_leaf_biomass_sum);
        External(INTERNODE_BIOMASS_SUM,&CulmStockModel::_internode_biomass_sum);
        External(PLANT_STOCK, &CulmStockModel::_plant_stock);
        External(PLANT_DEFICIT, &CulmStockModel::_plant_deficit);
        External(INTERNODE_DEMAND_SUM, &CulmStockModel::_internode_demand_sum);
        External(LEAF_DEMAND_SUM, &CulmStockModel::_leaf_demand_sum);
        External(INTERNODE_LAST_DEMAND_SUM,&CulmStockModel::_internode_last_demand_sum);
        External(LEAF_LAST_DEMAND_SUM, &CulmStockModel::_leaf_last_demand_sum);
        External(REALLOC_BIOMASS_SUM, &CulmStockModel::_realloc_biomass_sum);
        External(PLANT_PHASE, &CulmStockModel::_plant_phase);
        External(CULM_PHASE, &CulmStockModel::_culm_phase);
        External(PANICLE_DAY_DEMAND, &CulmStockModel::_panicle_day_demand);
        External(PANICLE_WEIGHT, &CulmStockModel::_panicle_weight);
        External(IS_FIRST_DAY_PI, &CulmStockModel::_is_first_day_pi);
        External(PEDUNCLE_DAY_DEMAND, &CulmStockModel::_peduncle_day_demand);
        External(PEDUNCLE_LAST_DEMAND, &CulmStockModel::_peduncle_last_demand);

    }

    virtual ~CulmStockModel()
    {}


    void compute(double t, bool /* update */) {
        if (_plant_phase == plant::INITIAL or _plant_phase == plant::VEGETATIVE) {
            _surplus = 0;
            return;
        }

        _stock_culm = _plant_stock * (_leaf_biomass_sum + _internode_biomass_sum) / _plant_biomass_sum;

        _deficit_culm =  _plant_deficit * (_leaf_biomass_sum + _internode_biomass_sum) / _plant_biomass_sum;

        _stock_internode = _maximum_reserve_in_internode * _internode_biomass_sum;

        _culm_biomass = _leaf_biomass_sum + _internode_biomass_sum;

        //MaxReservoirDispo
        _max_reservoir_dispo = (_maximum_reserve_in_internode *
                                _internode_biomass_sum) + (_leaf_stock_max * _leaf_biomass_sum);

        //CulmSupply
        _supply = _assim * _leaf_biomass_sum / _plant_leaf_biomass_sum;

        //Intermediate
        double stock = _plant_stock *
                (_leaf_biomass_sum + _internode_biomass_sum) /
                _plant_biomass_sum;
        double deficit = _plant_deficit *
                (_leaf_biomass_sum + _internode_biomass_sum) /
                _plant_biomass_sum;

        _intermediate = stock + deficit + _supply - _internode_demand_sum -
                _leaf_demand_sum - _leaf_last_demand_sum - _panicle_day_demand -
                _internode_last_demand_sum - _peduncle_day_demand - _peduncle_last_demand
                + _realloc_biomass_sum;

        //Deficit
        _deficit = std::min(0., _intermediate);

        //CulmSurplus
            if((_culm_phase != culm::INITIAL and _culm_phase != culm::VEGETATIVE) or (_plant_phase == plant::ELONG and _culm_phase == culm::VEGETATIVE)) {
                _surplus = std::max(0., _stock_culm - _internode_demand_sum -
                                    _leaf_demand_sum - _leaf_last_demand_sum -
                                    _internode_last_demand_sum - _panicle_day_demand
                                    - _peduncle_day_demand - _peduncle_last_demand
                                    + _supply - _max_reservoir_dispo +
                                    _realloc_biomass_sum);
            }

        if((_culm_phase != culm::INITIAL and _culm_phase != culm::VEGETATIVE) or (_plant_phase == plant::ELONG and _culm_phase == culm::VEGETATIVE)) {
            _stock = std::max(0., std::min(_max_reservoir_dispo, _intermediate));
        } else {
            _stock = _stock_culm;
        }

        _stock_internode = std::min(_stock, _maximum_reserve_in_internode * _internode_biomass_sum);
    }

    void init(double t, const ecomeristem::ModelParameters& parameters) {
        _parameters = parameters;
        //    parameters variables
        _maximum_reserve_in_internode = parameters.get("maximumReserveInInternode");
        _leaf_stock_max = parameters.get("leaf_stock_max");
        _realocationCoeff = parameters.get("realocationCoeff");

        //    computed variables (internal)
        _stock = 0;
        _supply = 0;
        _max_reservoir_dispo = 0;
        _intermediate = 0;
        _deficit = 0;
        _surplus = 0;
        _first_day = t;
        _leaf_stock = 0;
        _stock_culm = 0;
        _culm_biomass = 0;
        _stock_internode = 0;
        _deficit_culm = 0;
    }

private:
    ecomeristem::ModelParameters _parameters;

    // parameters
    double _maximum_reserve_in_internode;
    double _leaf_stock_max;
    double _realocationCoeff;

    //    internals - computed
    double _first_day;
    double _stock;
    double _supply;
    double _max_reservoir_dispo;
    double _intermediate;
    double _deficit;
    double _surplus;
    double _leaf_stock;
    double _stock_culm;
    double _culm_biomass;
    double _stock_internode;
    double _deficit_culm;

    //    externals
    plant::plant_phase _plant_phase;
    culm::culm_phase _culm_phase;
    double _assim;
    double _leaf_biomass_sum;
    double _internode_biomass_sum;
    double _plant_leaf_biomass_sum;
    double _plant_biomass_sum;
    double _plant_stock;
    double _plant_deficit;
    double _internode_demand_sum;
    double _leaf_demand_sum;
    double _internode_last_demand_sum;
    double _leaf_last_demand_sum;
    double _realloc_biomass_sum;
    double _panicle_day_demand;
    double _panicle_weight;
    double _last_leaf_biomass_sum;
    double _last_plant_biomass_sum;
    bool _is_first_day_pi;
    double _peduncle_day_demand;
    double _peduncle_last_demand;

};

} // namespace model

#endif // CULMSTOCKMODEL_H
