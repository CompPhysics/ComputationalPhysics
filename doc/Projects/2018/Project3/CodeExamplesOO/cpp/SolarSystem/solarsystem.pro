TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    vec3.cpp \
    celestialbody.cpp \
    system.cpp \
    newtoniangravity.cpp \
    forwardeuler.cpp

HEADERS += \
    vec3.h \
    celestialbody.h \
    system.h \
    newtoniangravity.h \
    forwardeuler.h
