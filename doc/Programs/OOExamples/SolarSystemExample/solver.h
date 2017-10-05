#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>
using std::vector;

class solver
{
public:
    friend class planet;

    // properties
    double radius,total_mass,G;
    int total_planets;
    vector<planet> all_planets;
    double totalKinetic;
    double totalPotential;

    // constants

    // initializers
    solver();
    solver(double radi);

    // functions
    void add(planet newplanet);
    void addM(planet newplanet);
    void GravitationalConstant();
    void print_position(std::ofstream &output, int dimension, double time, int number);
    void print_energy(std::ofstream &output, double time, double epsilon);
    void VelocityVerlet(int dimension, int integration_points, double final_time, int print_number, double epsilon);
    double **setup_matrix(int height, int width);
    void delete_matrix(double **matrix);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, double epsilon);
    void GravitationalForce_RK(double x_rel, double y_rel, double z_rel, double &Fx, double &Fy, double &Fz, double mass1, double mass2);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    double EnergyLoss();
    bool Bound(planet OnePlanet);

};

#endif // SOLVER_H
