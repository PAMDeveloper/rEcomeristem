QT       -= core gui
CONFIG += c++11
TEMPLATE = lib
CONFIG += static
NAME = rEcomeristem
SRC_ROOT = ../git/$$NAME/src
#R_PATH = C:/R
R_PATH = "C:/Program Files/R/R-3.3.0"
RCPP_PATH = $$R_PATH/library/Rcpp

CONFIG(debug, debug|release) {
    TARGET = $${NAME}d
} else {
    TARGET = $${NAME}
}

!contains(QT_ARCH, x86_64) {
    ARCHI = x86
    ARCHI_R = i386
} else {
    ARCHI = x64
    ARCHI_R = x64
}

CONFIG(static) {
    LINK = static
} else {
    LINK = shared
}

win32 {
    *-g++-64* {
        message("64bit here")
    }
    *-g++* {
        GCC_VERSION = $$system(gcc -dumpversion)
        equals(GCC_VERSION, "") {
            message("No compiler set: default 4.9.3")
            GCC_VERSION = 4.9.3
        }
        COMPILER = mingw-$$GCC_VERSION
    }
    *-msvc* {
        COMPILER = msvc14
    }
}

#macx {
#    *-g++* {}
#    *-xcode {}
#}

#unix{
#    *-g++* {}
#}


message($$TARGET - $$TEMPLATE - $$ARCHI - $$LINK - $$COMPILER)
DESTDIR = ../libs/$$COMPILER/$$ARCHI/$$LINK
message(to: $$DESTDIR)

LIBS += -L../ext_libs/$$COMPILER/$$ARCHI/static -llibpq \
        -L../libs/$$COMPILER/$$ARCHI/static -lartis -lecomeristem \
        -L$$RCPP_PATH/libs/$$ARCHI_R -lRcpp

INCLUDEPATH +=  ../ext_libs/include \
                $$SRC_ROOT \
                ../git/artis/src

INCLUDEPATH +=  ../../ext_libs/include \
                $$SRC_ROOT \
                ../git/artis/src \
                ../git/ecomeristem/src \
                $$RCPP_PATH/include \
                $$R_PATH/include

HEADERS += \
    $$SRC_ROOT/recomeristem_types.hpp

SOURCES += \
    $$SRC_ROOT/rcpp_ecomeristem.cpp \
    $$SRC_ROOT/RcppExports.cpp




