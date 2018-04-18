//#ifdef WITH_TRACE
#include "simpletrace.h"
SimpleTrace* SimpleTrace::lazy_singleton = nullptr;
SimpleTrace& SimpleTrace::trace() {
    if (lazy_singleton == nullptr) {
        lazy_singleton = new SimpleTrace();
    }
    return *lazy_singleton;
}
SimpleTrace::SimpleTrace(){}

std::map < std::string, int > KernelInfo::elt_dictionary;
std::vector < std::string > KernelInfo::elt_names;

//#endif
