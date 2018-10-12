#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include "vec3.h"


class CelestialBody
{
private:
public:
    // For now, public variables...
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;
    double kineticEnergy;
    double potentialEnergy;
    string name;

    // Default, parametrized constructor
    CelestialBody(vec3 newPosition, vec3 newVelocity, double newMass, string newName);

    // Destructor
    ~CelestialBody() {}


    void printObject();

    // No need for getters/setters when we have public variables!

};

#endif // CELESTIALBODY_H
