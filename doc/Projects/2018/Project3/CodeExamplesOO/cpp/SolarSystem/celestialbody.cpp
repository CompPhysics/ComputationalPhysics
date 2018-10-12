#include "celestialbody.h"

CelestialBody::CelestialBody()
{
}

CelestialBody::CelestialBody(vec3 newPosition, vec3 newVelocity, double newMass, string newName) :
    position(newPosition), velocity(newVelocity), mass(newMass), name(newName)
{
    /*
     * Writing the variable we take in as m_foo(foo) is the
     * same as setting it to m_foo = foo in the function body.
     */
    force = {0,0,0};
}

void CelestialBody::printObject()
{
    cout << name << " Position: " << position << " Velocity: " << velocity << endl;
}

