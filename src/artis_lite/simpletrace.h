#ifndef SIMPLETRACE_H
#define SIMPLETRACE_H

#include <vector>
#include <map>
#include <sstream>
using namespace std;

//#ifdef WITH_TRACE

namespace artis { namespace utils {
enum TraceType { NONE = 0, CHECK = 1, CONSTRUCT, SUBMODEL_ADD, INTERNAL_DECL,
                 EXTERNAL_DECL, INTERNAL_LINK, INIT, START, BEFORE_COMPUTE,
                 COMPUTE, PUT, AFTER_COMPUTE, DESTRUCT, KERNEL};
static const int DATE_FORMAT_YMD = 0;
}}

static const std::vector <std::string> TraceTypesStr = {
    "none", "check", "construct", "submodel_add", "internal_decl",
    "external_decl", "internal_link", "init", "start", "before_compute",
    "compute", "put", "after_compute", "destruct", "kernel"};


//class SimpleKernelInfo {
//public:
//    static int term(string s) {return 0;}
//    const std::string& var() const { return var; }
//    const std::string& value() const { return value; }
//    const std::string& tgt_model() const { return ""; }
//    const std::string& tgt_internal_var() const { return ""; }
//    unsigned int var_idx() const { return 0; }
//    unsigned int tgt_internal_var_idx() const { return 0; }
//    unsigned int tgt_model_idx() const { return 0; }
//    bool is_internal_var() const { return false; }
//    bool empty() const { return false; }
//    string to_string() const {return "";}
//};

class KernelInfo
{
public:
    static std::map < std::string, int > elt_dictionary;
    static std::vector < std::string > elt_names;
    static unsigned int term(const std::string& v) {
        if(elt_dictionary.find(v) == elt_dictionary.end()) {
            elt_dictionary[v] = elt_names.size();
            elt_names.push_back(v);
        }
        return elt_dictionary[v];
    }
    static const std::string & term(unsigned int i) {return elt_names[i];}


    KernelInfo() :
        _internal_var(false), _empty(true)
    {
        set_var("");
        set_value("");
        set_tgt_model("");
        set_tgt_internal_var("");
    }

    KernelInfo(std::string var, bool internal_var,
               std::string value) :
         _internal_var(internal_var), _empty(false)
    {
        set_var(var);
        set_value(value);
        set_tgt_model("");
        set_tgt_internal_var("");
    }

    KernelInfo(std::string var, const std::string tgt_model,
               std::string tgt_internal_var) :
        _internal_var(true), _empty(false)
    {
        set_var(var);
        set_value("");
        set_tgt_model(tgt_model);
        set_tgt_internal_var(tgt_internal_var);
    }

    KernelInfo(const std::string& var, bool internal_var) :
        _internal_var(internal_var), _empty(false)
    {
        set_var(var);
        set_value("");
        set_tgt_model("");
        set_tgt_internal_var("");
    }

    KernelInfo(std::string tgt_model) :
        _internal_var(false), _empty(false)
    {
        set_var("");
        set_value("");
        set_tgt_model(tgt_model);
        set_tgt_internal_var("");
    }

    void set_var(std::string v)
    { _var = KernelInfo::term(v); }

    void set_value(std::string v)
    { _value = KernelInfo::term(v); }

    void set_tgt_model(std::string v)
    { _tgt_model = KernelInfo::term(v); }

    void set_tgt_internal_var(std::string v)
    { _tgt_internal_var = KernelInfo::term(v); }


    const std::string& var() const
    { return KernelInfo::term(_var); }

    const std::string& value() const
    { return KernelInfo::term(_value); }

    const std::string& tgt_model() const
    { return KernelInfo::term(_tgt_model); }

    const std::string& tgt_internal_var() const
    { return KernelInfo::term(_tgt_internal_var); }

    unsigned int var_idx() const
    { return _var; }

    unsigned int tgt_internal_var_idx() const
    { return _tgt_internal_var; }

    unsigned int tgt_model_idx() const
    { return _tgt_model; }

    bool is_internal_var() const
    { return _internal_var; }

    bool empty() const
    { return _empty; }

    std::string to_string() const
    {
        std::ostringstream ss;
        if( !KernelInfo::term(_var).empty() ){
            ss << ":";
            if ( !_internal_var)
                ss << "*";
            ss << KernelInfo::term(_var);
        }
        if( !KernelInfo::term(_value).empty()) {
            ss << "=" << KernelInfo::term(_value);
        }
        if( !KernelInfo::term(_tgt_model).empty()) {
            ss << " -> " << KernelInfo::term(_tgt_model);
            if( !KernelInfo::term(_tgt_internal_var).empty()) {
                ss << ":" << KernelInfo::term(_tgt_internal_var);
            }
        }
        return ss.str();
    }

private:
    unsigned int _var;
    unsigned int _value;
    unsigned int _tgt_model;
    unsigned int _tgt_internal_var;
    bool        _internal_var;
    bool        _empty;
};


class SimpleTraceElement {
public:
    string _model_name; double _time; artis::utils::TraceType _type;
    KernelInfo _info;
    SimpleTraceElement(std::string model_name, double time, artis::utils::TraceType type, KernelInfo info) :
        _time(time), _type(type), _model_name(model_name), _info(info) {}

    double get_time() const { return _time; }
    const string& get_model_name() const { return _model_name; }
    artis::utils::TraceType get_type() const { return _type; }
    unsigned int get_model_name_idx() const {return 0;}
    const KernelInfo& get_kernel_info() const {return _info;}
    string to_string(int i = 0) const {
        std::ostringstream ss;
        ss << "TRACE " << _time << ": ";
        ss << "<" + TraceTypesStr[get_type()] + ">";
        ss << " " << get_model_name();
        if ( !get_kernel_info().empty()) {
            ss << get_kernel_info().to_string();
        }
        return ss.str();
    }
};

class SimpleTraceElements : public vector < SimpleTraceElement > {
public:
    string to_string(int i = 0) const {
        std::ostringstream ss;
        for (auto it = SimpleTraceElements::begin(); it != SimpleTraceElements::end(); ++it)
            ss << it->to_string(i) << std::endl;
        return ss.str();
    }
};

class SimpleTrace {
public:
    static SimpleTrace& trace();
    void clear(){_trace.clear();}
    const SimpleTraceElements & elements() {return _trace;}
    void addElement(SimpleTraceElement t) {_trace.push_back(t);}
private:
    static SimpleTrace* lazy_singleton;
    SimpleTrace();
    SimpleTraceElements _trace;
};

//#endif

#endif // SIMPLETRACE_H
