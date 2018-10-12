#include <iostream>
#include <cmath>
#include "vec3.h"
#include "celestialbody.h"
#include "system.h"
#include "forwardeuler.h"
#include "newtoniangravity.h"

using std::cout;
using std::endl;

int main()
{
    // Some basic operations...
    vec3 v1(1,1,1);
    vec3 v2(2,2,2);
    v1.print();
    v2.print();
    vec3 v3, v4;
    v3 = v1;
    v3.print();
    v4 = v1 + v2 - v3;
    v4.print();


    // Project begins here
    // Initializes sun, earth, jupiter
    double massSun = 1988500e24; // kg
    vec3 posSun(0,0,0);
    vec3 velSun(0,0,0);

    double massEarth = 5.97219e24; // massEarth / massSun
    vec3 posEarth(1,0,0);
    vec3 velEarth(0,2*M_PI,0);

    double massJupiter = 1.0;
    vec3 posJupiter(0,0,0);
    vec3 velJupiter(0,0,0);

    CelestialBody *sun = new CelestialBody(posSun, velSun, massSun, "sun");
    CelestialBody *earth = new CelestialBody(posEarth, velEarth, massEarth, "earth");
    CelestialBody *jupiter = new CelestialBody(posJupiter, velJupiter, massJupiter, "jupiter");

    // Specifies the number of steps, time to run for and step size
    unsigned long NSteps = 100;
    double T = 1; // Time, years
    double h = T / double(NSteps);
    double G = 4*M_PI*M_PI;

    System S;
    ForwardEuler integrator;
    NewtonianGravity force(G, massSun);

    S.addObject(sun);
    S.addObject(earth);
//    S.addObject(jupiter); // First verify that the 2-body problem produces the correct results

    cout << "Number of bodies in the solar system is: " << S.bodies.size() << endl;

    // Prints start position
    S.bodies[1]->printObject();

    for (int iStep = 0; iStep < NSteps; iStep++)
    {
        S.resetForces(); // Since we only stores one step at the time...
        force.calculateForces(&S);
        integrator.integrate(&S, h);

        for (int iObj = 0; iObj < S.bodies.size(); iObj++)
        {
            S.bodies[iObj]->printObject();
        }

    }

    // Prints end position, which should be approximately the same as the start
    S.bodies[1]->printObject();


    return 0;
}
