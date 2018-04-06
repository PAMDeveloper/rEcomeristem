#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H

#include <vector>
#include <iterator>
#include <list>
#include <ModelParameters.hpp>
using namespace std;

//useless - for signature compatibility
namespace artis { namespace kernel {
    enum ValueTypeID { DOUBLE, INT, BOOL, DOUBLE_VECTOR, INT_VECTOR, BOOL_VECTOR, STRING, USER };
}}


template < typename T >
class SimpleView {
    typedef vector < unsigned int > Selector;

public:
    void selector(const string& name,
                  artis::kernel::ValueTypeID value_type,
                  const Selector& chain)
    {}
};

class SimpleContext {
public:
    double start, end;
    SimpleContext(double start, double end):
        start(start), end(end) {}
};

class AbstractSimpleModel {
public:
    virtual void compute(double t, bool update = false) = 0;
    virtual void init(double t, const ecomeristem::ModelParameters& parameters) = 0;
    virtual void operator()(double t) { this->compute(t); }
};

template < typename T, typename U, typename V >
class SimpleSimulator {
public:
    T * _model;
    SimpleSimulator(T * model, U parameters) : _model(model) {}

    void init(double time, const V& parameters)
    { _model->init(time, parameters); }

    void run(const SimpleContext & context)
    {}

    void attachView(const string& name, SimpleView < V > * view)
    { }

};


template < typename T >
class SimpleModel : public AbstractSimpleModel {

protected:
    double last_time;
public:
    vector<double T::*> i_double;
    vector<int T::*> i_int;
    vector<bool T::*> i_bool;
    vector<plant::plant_state T::*> i_plant_state;
    vector<plant::plant_phase T::*> i_plant_phase;
    vector<culm::culm_phase T::*> i_culm_phase;
    vector<internode::internode_phase T::*> i_internode_phase;
    vector<leaf::leaf_phase T::*> i_leaf_phase;
    vector<peduncle::peduncle_phase T::*> i_peduncle_phase;

    vector<double T::*> e_double;
    vector<int T::*> e_int;
    vector<bool T::*> e_bool;
    vector<plant::plant_state T::*> e_plant_state;
    vector<plant::plant_phase T::*> e_plant_phase;
    vector<culm::culm_phase T::*> e_culm_phase;
    vector<internode::internode_phase T::*> e_internode_phase;
    vector<leaf::leaf_phase T::*> e_leaf_phase;
    vector<peduncle::peduncle_phase T::*> e_peduncle_phase;

    vector< AbstractSimpleModel * > subModels;
//    ecomeristem::ModelParameters * parameters;

    template <typename W>
    void add_vec(unsigned int index, W T::* var, vector<W T::*> & vec) {
        if(vec.size() <= index)
            while(vec.size() <= index)
                vec.push_back(var);
        else vec[index] = var;
    }

    void internal_(unsigned int index, const string& /*n*/, double T::* var) {add_vec<double>(index,var,i_double);}
    void internal_(unsigned int index, const string& /*n*/, int T::* var) {add_vec<int>(index,var,i_int);}
    void internal_(unsigned int index, const string& /*n*/, bool T::* var) {add_vec<bool>(index,var,i_bool);}
    void internal_(unsigned int index, const string& /*n*/, plant::plant_state T::* var) {add_vec<plant::plant_state>(index,var,i_plant_state);}
    void internal_(unsigned int index, const string& /*n*/, plant::plant_phase T::* var) {add_vec<plant::plant_phase>(index,var,i_plant_phase);}
    void internal_(unsigned int index, const string& /*n*/, culm::culm_phase T::* var) {add_vec<culm::culm_phase>(index,var,i_culm_phase);}
    void internal_(unsigned int index, const string& /*n*/, internode::internode_phase T::* var) {add_vec<internode::internode_phase>(index,var,i_internode_phase);}
    void internal_(unsigned int index, const string& /*n*/, leaf::leaf_phase T::* var) {add_vec<leaf::leaf_phase>(index,var,i_leaf_phase);}
    void internal_(unsigned int index, const string& /*n*/, peduncle::peduncle_phase T::* var) {add_vec<peduncle::peduncle_phase>(index,var,i_peduncle_phase);}

    void external_(unsigned int index, const string& /*n*/, double T::* var) {add_vec<double>(index,var,e_double);}
    void external_(unsigned int index, const string& /*n*/, int T::* var) {add_vec<int>(index,var,e_int);}
    void external_(unsigned int index, const string& /*n*/, bool T::* var) {add_vec<bool>(index,var,e_bool);}
    void external_(unsigned int index, const string& /*n*/, plant::plant_state T::* var) {add_vec<plant::plant_state>(index,var,e_plant_state);}
    void external_(unsigned int index, const string& /*n*/, plant::plant_phase T::* var) {add_vec<plant::plant_phase>(index,var,e_plant_phase);}
    void external_(unsigned int index, const string& /*n*/, culm::culm_phase T::* var) {add_vec<culm::culm_phase>(index,var,e_culm_phase);}
    void external_(unsigned int index, const string& /*n*/, internode::internode_phase T::* var) {add_vec<internode::internode_phase>(index,var,e_internode_phase);}
    void external_(unsigned int index, const string& /*n*/, leaf::leaf_phase T::* var) {add_vec<leaf::leaf_phase>(index,var,e_leaf_phase);}
    void external_(unsigned int index, const string& /*n*/, peduncle::peduncle_phase T::* var) {add_vec<peduncle::peduncle_phase>(index,var,e_peduncle_phase);}

    template < typename W > W get(double /*t*/, unsigned int index) { return W(); }
    template < typename W, typename U > W get(double t, unsigned int index) { return get<W>(t,index);}
    template <> double get<double>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<double SimpleModel<T>::*>(i_double[index]);}
    template <> int get<int>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<int SimpleModel<T>::*>(i_int[index]);}
    template <> bool get<bool>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<bool SimpleModel<T>::*>(i_bool[index]);}
    template <> plant::plant_state get<plant::plant_state>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<plant::plant_state SimpleModel<T>::*>(i_plant_state[index]);}
    template <> plant::plant_phase get<plant::plant_phase>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<plant::plant_phase SimpleModel<T>::*>(i_plant_phase[index]);}
    template <> culm::culm_phase get<culm::culm_phase>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<culm::culm_phase SimpleModel<T>::*>(i_culm_phase[index]);}
    template <> internode::internode_phase get<internode::internode_phase>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<internode::internode_phase SimpleModel<T>::*>(i_internode_phase[index]);}
    template <> leaf::leaf_phase get<leaf::leaf_phase>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<leaf::leaf_phase SimpleModel<T>::*>(i_leaf_phase[index]);}
    template <> peduncle::peduncle_phase get<peduncle::peduncle_phase>(double /*t*/, unsigned int index)
    {return static_cast<SimpleModel<T>*>(this)->*static_cast<peduncle::peduncle_phase SimpleModel<T>::*>(i_peduncle_phase[index]);}

    template < typename W > void put(double /*t*/, unsigned int index, W value) {}
    template <> void put<double>(double /*t*/, unsigned int index, double value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<double SimpleModel<T>::*>(e_double[index]) = value;}
    template <> void put<int>(double /*t*/, unsigned int index, int value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<int SimpleModel<T>::*>(e_int[index]) = value;}
    template <> void put<bool>(double /*t*/, unsigned int index, bool value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<bool SimpleModel<T>::*>(e_bool[index]) = value;}
    template <> void put<plant::plant_state>(double /*t*/, unsigned int index, plant::plant_state value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<plant::plant_state SimpleModel<T>::*>(e_plant_state[index]) = value;}
    template <> void put<plant::plant_phase>(double /*t*/, unsigned int index, plant::plant_phase value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<plant::plant_phase SimpleModel<T>::*>(e_plant_phase[index]) = value;}
    template <> void put<culm::culm_phase>(double /*t*/, unsigned int index, culm::culm_phase value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<culm::culm_phase SimpleModel<T>::*>(e_culm_phase[index]) = value;}
    template <> void put<internode::internode_phase>(double /*t*/, unsigned int index, internode::internode_phase value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<internode::internode_phase SimpleModel<T>::*>(e_internode_phase[index]) = value;}
    template <> void put<leaf::leaf_phase>(double /*t*/, unsigned int index, leaf::leaf_phase value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<leaf::leaf_phase SimpleModel<T>::*>(e_leaf_phase[index]) = value;}
    template <> void put<peduncle::peduncle_phase>(double /*t*/, unsigned int index, peduncle::peduncle_phase value)
    {static_cast<SimpleModel<T>*>(this)->*static_cast<peduncle::peduncle_phase SimpleModel<T>::*>(e_peduncle_phase[index]) = value;}

    //voir les submodel internals
    void link_internal_(unsigned int index, string var_name, AbstractSimpleModel * model, int sub_index, string sub_var_name)
    {}

    void setsubmodel(unsigned int index, AbstractSimpleModel * model)
    {}

};

#define DOUBLEESCAPE(a) #a
#define ESCAPEQUOTE(a) DOUBLEESCAPE(a)
#define Internal(index, var) internal_(index, string(ESCAPEQUOTE(index)), var)
#define External(index, var) external_(index, string(ESCAPEQUOTE(index)), var)
#define InternalS(index, var, sub_index) link_internal_(index, string(ESCAPEQUOTE(index)), var, sub_index, string(ESCAPEQUOTE(index)))


#endif // SIMPLEMODEL_H
