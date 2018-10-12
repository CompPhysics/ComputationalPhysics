#include "system.h"

System::System()
{

}

void System::resetForces() {
    /*
     * Ensures the force vector in each object is set to zero.
     */
    for (CelestialBody *obj : bodies) {
        obj->force = {0,0,0};
    }
}
