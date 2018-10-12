#ifndef FORWARDEULER_H
#define FORWARDEULER_H

#include "system.h"

class ForwardEuler
{
public:
    ForwardEuler();
    ~ForwardEuler() {}

    void integrate(System *system, const double h);
};

#endif // FORWARDEULER_H
