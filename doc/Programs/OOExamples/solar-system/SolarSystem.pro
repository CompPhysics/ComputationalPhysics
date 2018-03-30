TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    solarsystem.cpp \
    vec3.cpp \
    euler.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h \
    vec3.h \
    euler.h

