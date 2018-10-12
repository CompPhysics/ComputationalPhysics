#ifndef NEWTONIANGRAVITY_H
#define NEWTONIANGRAVITY_H

#include <vector>
#include "system.h"

using std::vector;

class NewtonianGravity
{
private:
    double m_G;
    double m_sunMass;
public:
    NewtonianGravity(double G, double sunMass);
    ~NewtonianGravity() {}

    void calculateForces(System *system);
};

#endif // NEWTONIANGRAVITY_H
