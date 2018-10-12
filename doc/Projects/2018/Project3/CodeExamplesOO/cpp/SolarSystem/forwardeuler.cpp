#include "forwardeuler.h"

ForwardEuler::ForwardEuler()
{

}

void ForwardEuler::integrate(System *system, const double h)
{
    /*
     * A basic forward Euler(or Euler Cromer to be precise)
     */

    for (CelestialBody *obj : system->bodies) {

        obj->velocity += obj->force*h / obj->mass;

        obj->position += obj->velocity*h;
    }
}
