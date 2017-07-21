#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <artis/kernel/AbstractAtomicModel.hpp>
#include <artis/kernel/AbstractCoupledModel.hpp>
#include <artis/kernel/Simulator.hpp>
#include <artis/observer/Observer.hpp>
#include <artis/observer/View.hpp>
#include <artis/utils/DoubleTime.hpp>
#include <ModelParameters.hpp>
#include <utils/julianconverter.h>

#include <type_traits> //necessary for enum type checking on states

class PlantModel;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace peduncle {
enum peduncle_phase {   INITIAL = 0,
                        TRANSITION = 1,
                        REALIZATION = 2,
                        FLO = 3,
                        END_FILLING = 4,
                        MATURITY = 5,
                        DEAD = 6 };

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

enum plant_phase {  INITIAL = 0,
                    VEGETATIVE = 1,
                    ELONG = 2,
                    PI = 3,
                    PRE_FLO = 4,
                    FLO = 5,
                    END_FILLING = 6,
                    MATURITY = 7,
                    DEAD = 8 };

inline bool operator&(plant_state a, plant_state b)
{return (static_cast<int>(a) & static_cast<int>(b)) != 0;}
inline void operator<<(plant_state& a, plant_state b)
{a = static_cast<plant_state>(static_cast<int>(a) | static_cast<int>(b));}
inline void operator>>(plant_state& a, plant_state b)
{a = static_cast<plant_state>(static_cast<int>(a) & !static_cast<int>(b));}

//enum plant_phase { INIT = 0,
//                   INITIAL = 1,
//                   GROWTH = 2,
//                   NOGROWTH = 3,
//                   NEW_PHYTOMER = 5,
//                   NOGROWTH2 = 18,
//                   NOGROWTH3 = 19,
//                   NOGROWTH4 = 20,
//                   NEW_PHYTOMER3 = 23,
//                   LIG = 24,
//                   KILL = 25 };

//enum plant_state {  VEGETATIVE = 0,
//                    PRE_ELONG,
//                    ELONG,
//                    PRE_PI,
//                    PI,
//                    PRE_FLO,
//                    FLO,
//                    END_FILLING,
//                    MATURITY,
//                    DEAD };
}


struct GlobalParameters
{ };

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

#endif // DEFINES_HPP
