#include "newtoniangravity.h"
#include <cmath>

NewtonianGravity::NewtonianGravity(double G, double sunMass)
{
    m_G = G;
    m_sunMass = sunMass;
}


void NewtonianGravity::calculateForces(System *system)
{
    /*
     * Takes the system class and calculates forces between the objects
     */


}
