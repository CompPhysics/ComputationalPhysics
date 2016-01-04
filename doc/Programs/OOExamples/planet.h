#ifndef PLANET_H
#define PLANET_H

class planet
{
public:

    double position[3];
    double velocity[3];
    double mass;

    planet(double mas, double x,double y, double z, double vx, double vy, double vz);
    planet();

};

#endif // PLANET_H
