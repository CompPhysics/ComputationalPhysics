#include "system.h"

System::System()
{

}

void System::resetForces() {
    /*
     * Ensures the force vector in each object is set to zero.
     */

    // This is the same as for element in list in python!
    // One of the many handy methods in C++11.
    for (CelestialBody *obj : bodies) {
        obj->force = {0,0,0};
    }
    // It is the same as following:
//    for (int iObj = 0; iObj < bodies.size(); iObj++) {
//        bodies[iObj]->force = {0,0,0};
//    }

    // Note that vec3 allows setting vectors with my_vec = {x,y,z}
}
