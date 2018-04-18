#ifndef DEFINES_HPP
#define DEFINES_HPP

#include "iso646.h"
#include <type_traits> //necessary for enum type checking on states

class PlantModel;


struct GlobalParameters
{ };


#ifdef UNSAFE_RUN

namespace peduncle {
typedef int peduncle_phase;
static const unsigned int INITIAL = 0;
static const unsigned int TRANSITION = 1;
static const unsigned int REALIZATION = 2;
static const unsigned int FLO = 3;
static const unsigned int END_FILLING = 4;
static const unsigned int MATURITY = 5;
static const unsigned int DEAD = 6;
}

namespace internode {
typedef int internode_phase;
static const unsigned int INITIAL = 0;
static const unsigned int VEGETATIVE = 1;
static const unsigned int REALIZATION = 2;
static const unsigned int MATURITY = 3;
static const unsigned int DEAD = 4;
}

namespace leaf {
typedef int leaf_phase;
static const unsigned int INITIAL = 0;
static const unsigned int VEGETATIVE = 1;
static const unsigned int LIG = 2;
static const unsigned int DEAD = 3;
}

namespace culm {
typedef int culm_phase;
static const unsigned int INITIAL = 0;
static const unsigned int VEGETATIVE = 1;
static const unsigned int ELONG = 2;
static const unsigned int PRE_PI = 3;
static const unsigned int PI = 4;
static const unsigned int PRE_FLO = 5;
static const unsigned int FLO = 6;
static const unsigned int END_FILLING = 7;
static const unsigned int MATURITY = 8;
static const unsigned int DEAD = 9;
}

struct flag {
    flag() = default;
    flag(int v): val(v){}
    int val;
    inline flag& operator=(const flag& b) {val = b.val; return *this;}
    inline flag& operator=(const unsigned int& b) {val = b; return *this;}
    inline operator int() const {return val;}
    inline explicit operator int*() const { return nullptr; }
};
inline bool operator&(const flag& a, const unsigned int& b) {return (a.val & b) != 0;}
inline void operator<<(flag& a, const unsigned int& b) {a.val = a.val | b;}
inline void operator>>(flag& a, const unsigned int& b) {a.val = a.val & !b;}

// Plant enums
namespace plant {
typedef int plant_phase;
static const unsigned int INITIAL = 0;
static const unsigned int VEGETATIVE = 1;
static const unsigned int ELONG = 2;
static const unsigned int PI = 3;
static const unsigned int PRE_FLO = 4;
static const unsigned int FLO = 5;
static const unsigned int END_FILLING = 6;
static const unsigned int MATURITY = 7;
static const unsigned int DEAD = 8;

typedef flag plant_state;
static const unsigned int NO_STATE = 0;
static const unsigned int NOGROWTH = 1;
static const unsigned int NEW_PHYTOMER_AVAILABLE = 2;
static const unsigned int LIG = 4;
static const unsigned int KILL = 8;

}

#include <utils/juliancalculator.h>
#include <artis_lite/simplemodel.h>
#include <memory>
#include <deque>
#include <limits>

namespace artis { namespace utils { namespace DateTime {
static string toJulianDayFmt(double, int) {
    return "";
}
}}}

struct DoubleTime {
    static constexpr double negative_infinity = -numeric_limits < double >::infinity();
    static constexpr double positive_infinity = numeric_limits < double >::infinity();
    static constexpr double null = 0;
};


using Model = SimpleModel < ecomeristem::ModelParameters >;
template < typename T > using AtomicModel = SimpleModel < T >;
template < typename T > using CoupledModel = SimpleModel < T >;

typedef SimpleSimulator < PlantModel, GlobalParameters, ecomeristem::ModelParameters > EcomeristemSimulator;
typedef SimpleContext EcomeristemContext;
typedef SimpleView View;
typedef SimpleObserver Observer;

//#ifdef WITH_TRACE
using Trace = SimpleTrace;
//using KernelInfo = SimpleKernelInfo;
template < class T > using TraceElement = SimpleTraceElement;
template < class T > using TraceElements = std::vector < TraceElement <T> >;
//#endif

#else

#include <utils/juliancalculator.h>
#include <ModelParameters.hpp>
#include <artis/kernel/Simulator.hpp>
#include <artis/kernel/AbstractAtomicModel.hpp>
#include <artis/utils/Trace.hpp>

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
enum internode_phase {  INITIAL = 0,
                        VEGETATIVE = 1,
                        REALIZATION = 2,
                        MATURITY = 3,
                        DEAD = 4 };
}

namespace leaf {
enum leaf_phase   { INITIAL = 0,
                    VEGETATIVE = 1,
                    LIG = 2,
                    DEAD = 3 };
}

namespace culm {
enum culm_phase {   INITIAL = 0,
                    VEGETATIVE = 1,
                    ELONG = 2,
                    PRE_PI = 3,
                    PI = 4,
                    PRE_FLO = 5,
                    FLO = 6,
                    END_FILLING = 7,
                    MATURITY = 8,
                    DEAD = 9 };

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
