#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <ModelParameters.hpp>
#include <utils/julianconverter.h>
#include <type_traits> //necessary for enum type checking on states

class PlantModel;

//#ifndef M_PI
//#define M_PI 3.14159265358979323846
//#endif

namespace peduncle {
enum peduncle_phase {   INITIAL = 0,
                        TRANSITION = 1,
                        REALIZATION = 2,
                        FLO = 3,
                        END_FILLING = 4,
                        MATURITY = 5,
                        DEAD = 6 };

}

namespace internode {
enum internode_phase {  INITIAL, VEGETATIVE, REALIZATION,
                        MATURITY, DEAD };
}

namespace leaf {
enum leaf_phase   { INITIAL, VEGETATIVE, LIG, DEAD };
}

namespace culm {
enum culm_phase {   INITIAL = 0,
                    VEGETATIVE = 1,
                    ELONG = 3,
                    PRE_PI = 4,
                    PI = 5,
                    PRE_FLO = 6,
                    FLO = 7,
                    END_FILLING = 8,
                    MATURITY = 9,
                    DEAD = 10 };

}

// Plant enums
namespace plant {
enum plant_state { NO_STATE = 0,
                   NOGROWTH = 1,
                   NEW_PHYTOMER_AVAILABLE = 2,
                   LIG = 4,
                   KILL = 8 };

enum plant_phase {  INITIAL = 0,
                    VEGETATIVE = 1,
                    ELONG = 2,
                    PI = 3,
                    PRE_FLO = 4,
                    FLO = 5,
                    END_FILLING = 6,
                    MATURITY = 7,
                    DEAD = 8 };

template <typename E, typename std::enable_if<std::is_enum<E>::value>::type* = nullptr>
class States
{
public:
    int states;

    bool operator&(E state)
    {return (states & static_cast<int>(state)) != 0;}
    void operator<<(E state)
    {states = states | static_cast<int>(state);}
    void operator>>(E state)
    {states = states & !static_cast<int>(state);}

    operator int() {return states;}
    bool is(E state)
    {return (states & static_cast<int>(state)) != 0;}
    void add(E state)
    {states = states | static_cast<int>(state);}
    void del(E state)
    {states = states & !static_cast<int>(state);}
};

typedef States <plant_state> plant_states;

inline bool operator&(plant_state a, plant_state b)
{return (static_cast<int>(a) & static_cast<int>(b)) != 0;}
inline void operator<<(plant_state& a, plant_state b)
{a = static_cast<plant_state>(static_cast<int>(a) | static_cast<int>(b));}
inline void operator>>(plant_state& a, plant_state b)
{a = static_cast<plant_state>(static_cast<int>(a) & !static_cast<int>(b));}
}


struct GlobalParameters
{ };


#if 1
#include <artis_lite/simplemodel.h>
#include <memory>
#include <deque>


using Model = SimpleModel < ecomeristem::ModelParameters >;

template < typename T >
using AtomicModel = SimpleModel < T >;

template < typename T >
using CoupledModel = SimpleModel < T >;

typedef SimpleSimulator < PlantModel, GlobalParameters, ecomeristem::ModelParameters > EcomeristemSimulator;

typedef SimpleContext EcomeristemContext;

typedef SimpleView < ecomeristem::ModelParameters > View;

//typedef artis::observer::Observer < artis::utils::DoubleTime,
//                                    ecomeristem::ModelParameters > Observer;

//using Trace = artis::utils::Trace < artis::utils::DoubleTime >;

//using TraceElement = artis::utils::TraceElement < artis::utils::DoubleTime >;


#else

#include <artis/kernel/AbstractAtomicModel.hpp>
#include <artis/kernel/AbstractCoupledModel.hpp>
#include <artis/kernel/Simulator.hpp>
#include <artis/observer/Observer.hpp>
#include <artis/observer/View.hpp>
#include <artis/utils/DoubleTime.hpp>

using Model = artis::kernel::AbstractModel < artis::utils::DoubleTime,
                                             ecomeristem::ModelParameters >;

using Trace = artis::utils::Trace < artis::utils::DoubleTime >;

using TraceElement = artis::utils::TraceElement < artis::utils::DoubleTime >;

template < typename T >
using AtomicModel = artis::kernel::AbstractAtomicModel <
    T, artis::utils::DoubleTime, ecomeristem::ModelParameters >;

template < typename T >
using CoupledModel = artis::kernel::AbstractCoupledModel <
    T, artis::utils::DoubleTime, ecomeristem::ModelParameters, GlobalParameters >;

typedef artis::observer::Observer < artis::utils::DoubleTime,
                                    ecomeristem::ModelParameters > Observer;

typedef artis::observer::View < artis::utils::DoubleTime,
                                ecomeristem::ModelParameters > View;

typedef artis::kernel::Simulator < PlantModel,
                                   artis::utils::DoubleTime,
                                   ecomeristem::ModelParameters,
                                   GlobalParameters > EcomeristemSimulator;

typedef artis::context::Context < artis::utils::DoubleTime > EcomeristemContext;
#endif

#endif // DEFINES_HPP
