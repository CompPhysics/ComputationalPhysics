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

    // Default constructor
    CelestialBody();

    // Destructor
    ~CelestialBody() {}

    // Default, parametrized constructor
    CelestialBody(vec3 newPosition, vec3 newVelocity, double newMass, string newName);

    void printObject();

    // No need for getters/setters when we have public variables!
//    // Getters
//    string getName() { return m_name; }
//    vec3 getPosition() { return m_position; }
//    vec3 getVelocity() { return m_velocity; }
//    double getMass() { return m_mass; }
//    double getKineticEnergy() { return m_kineticEnergy; }
//    double getPotentialEnergy() { return m_potentialEnergy; }

//    // Setters
//    void setName(string name) { m_name = name; }
//    void setPosition(vec3 position) { m_position = position; }
//    void setVelocity(vec3 velocity) { m_velocity = velocity; }
//    void setMass(double mass) { m_mass = mass; }
//    void setKineticEnergy(double kineticEnergy) { m_kineticEnergy = kineticEnergy; }
//    void setPotentialEnergy(double potentialEnergy) { m_potentialEnergy = potentialEnergy; }

};

#endif // CELESTIALBODY_H
