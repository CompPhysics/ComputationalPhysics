#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include "planet.h"
#include "solver.h"

using namespace std;

void print_initial(int dimension, double time_step, double final_time, double *x_initial, double *v_initial, int N);
void print_final(int dimension, double *x_final, double *v_final);
void randomUniformSphere(double R0,double &x,double &y,double &z, default_random_engine *generator);
void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator);

int main()
{    
    // WHICH SYSTEM WOULD YOU LIKE TO RUN?

    // 1) Analytical test with box on a spring; set spring_test = true
    bool spring_test = false;

    // 2) Two planets in a gravitational field; set binary = true
    bool binary = true;

    // 3) Compare RK4 and VV in a full planet cluster, set RK4vsVV = true
    bool RK4vsVV = false;

    // 4) Full planet cluster in a gravitational field for ONE numerical method; set cluster = true
    bool cluster = false;
    bool constant = false; // The system has constant mass


    // Numerical setup
    int integration_points;  // No. of integration points
    double final_time;       // End time of calculation
    int dimension;           // No. of spatial dimensions

    bool force;             // false = run program with analytical spring force, true = run program with gravitational potential
    bool simple;            // true = Sun-Earth-like system with G = 4*pi*pi instead of dimensionless G used for cluster later


    // 1) ANALYTIC: Spring force
    if(spring_test){
        cout << "ANALYTICAL" << endl;
        dimension = 1;

        force = false; // false = run program with analytical spring force
        simple = true; // true = Sun-Earth-like system with G = 4*pi*pi

        integration_points = 100;
        final_time = 100;
        cout << "Time step: " << final_time/((double) integration_points) << ", integration points: " << integration_points << endl;

        // RK4 test
        planet planet1(1.,1.,0,0,0,0,0); // (mass,x,y,z,vx,vy,vz)
        solver testRK;
        testRK.add(planet1);
        testRK.RungeKutta4(dimension,integration_points,final_time,force,simple,1,0.);

        // VV test
        planet planet2(1.,1.,0,0,0,0,0); // (mass,x,y,z,vx,vy,vz)
        solver testVV;
        testVV.add(planet2);
        testVV.VelocityVerlet(dimension,integration_points,final_time,force,simple,1,0.);
    }

    // 2) Binary planets
    if(binary){
        cout << "BINARY SYSTEM" << endl;
        dimension = 3;
        force = true; // true = run program with gravitational potential
        simple = true; // true = Sun-Earth-like system with G = 4*pi*pi

        integration_points = 10000;
        final_time = 50.;

        double time_step = final_time/((double) integration_points);
        double x[3],v[3];

        planet planet1(0.000003,1.,0.0,0.0,0.0,6.3,0.); // Earth: (mass,x,y,z,vx,vy,vz)
        planet planet2(1.,0.,0.,0.,0.,0.,0.);           // Sun: (mass,x,y,z,vx,vy,vz)

        // RK4
        solver binary_rk(5.0);
        binary_rk.add(planet1);
        binary_rk.add(planet2);

        for(int j=0;j<dimension;j++){
            x[j] = planet1.position[j];
            v[j] = planet1.velocity[j];
        }

        // VV
        solver binary_vv(5.0);
        binary_vv.add(planet1);
        binary_vv.add(planet2);

        print_initial(dimension,time_step,final_time,x,v,integration_points);

        // Evolution of binary system
        cout << endl << "RK4: " << endl;
        binary_rk.RungeKutta4(dimension,integration_points,final_time,force,simple,1,0.);

        for(int j=0;j<dimension;j++){
            x[j] = binary_rk.all_planets[0].position[j];
            v[j] = binary_rk.all_planets[0].velocity[j];
        }
        print_final(dimension,x,v);

        cout << endl << "VV:" << endl;
        binary_vv.VelocityVerlet(dimension,integration_points,final_time,force,simple,1,0.);

        for(int j=0;j<dimension;j++){
            x[j] = binary_vv.all_planets[0].position[j];
            v[j] = binary_vv.all_planets[0].velocity[j];
        }
        print_final(dimension,x,v);
    }

    // 3) N-BODY CASE: COMPARE METHODS
    if(RK4vsVV){
        cout << "CLUSTER" << endl;
        dimension = 3;
        force = true;           // true = run program with gravitational potential
        simple = false;         // false = use dimensionless G

        integration_points = 1000;
        final_time = 1.;        // in units of t_crunch
        double R0 = 20.;        // Radius of solver, in units of lightyears
        int objects = 100;      // Number of planets to be added in solver

        cout << "Time step: " << final_time/((double) integration_points) << endl;
        cout << "Integration points: " << integration_points << endl;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);

        // initialize mass and position, to be randomly distributed
        double m,x,y,z;
        m = x = y = z = 0.0;

        // mean and standard deviation of stellar mass
        double mean = 10.;
        double deviation = 1.;

        solver MM15_rk(R0); // set up cluster for RK4
        solver MM15_vv(R0); // set up cluster for VV

        // add planets to the two systems
        for(int i=0;i<objects;i++){
            Gaussian_distribution(mean,deviation,m,&generator);
            randomUniformSphere(R0,x,y,z,&generator);
            planet planeti(m,x,y,z,0,0,0);
            MM15_rk.add(planeti);
            MM15_vv.add(planeti);
        }
        cout << "The planet cluster MM15 contains " << MM15_rk.total_planets << " planet(s)." << endl;

        // run system through RK4/VV, all data is written to file as we go
        cout << "RK4" << endl;
        MM15_rk.RungeKutta4(dimension,integration_points,final_time,force,simple,objects,0.);
        cout << "VV" << endl;
        MM15_vv.VelocityVerlet(dimension,integration_points,final_time,force,simple,objects,0.);
    }

    // 4) SOLVER (PLANET CLUSTER) MODEL
    if(cluster){
        cout << "CLUSTER" << endl;
        dimension = 3;
        force = true;   // true = run program with gravitational potential
        simple = false; // false = use dimensionless G

        integration_points = 10000;
        final_time = 20.; // in units of t_crunch

        double epsilon = 0.1; // Smoothing factor

        cout << "Time step: " << final_time/((double) integration_points) << endl;
        cout << "Integration points: " << integration_points << endl;
        cout << "epsilon = " << epsilon << endl;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);

        double R0 = 20.;        // Radius of solver, in units of lightyears
        int objects = 200;      // Number of planets to be added in solver

        // initialize mass and position, to be randomly distributed
        double m,x,y,z;
        m = x = y = z = 0.0;

        // mean and standard deviation of stellar mass
        double mean = 10.;
        double deviation = 1.;

        solver MM15(R0); // set up cluster with radius R0

        // add planets to the two systems
        for(int i=0;i<objects;i++){
            Gaussian_distribution(mean,deviation,m,&generator);
            randomUniformSphere(R0,x,y,z,&generator);
            planet planeti(m,x,y,z,0,0,0);
            MM15.add(planeti);
        }
        cout << "The planet cluster MM15 contains " << MM15.total_planets << " planet(s)." << endl;

        // if-test to keep the mass of the system constant for changing values of N
        if(constant){
            double const_mass = 1000.;
            double total = MM15.total_mass;
            double current_mass = 0.;
            for(int nr=0;nr<MM15.total_planets;nr++){
                planet &Current = MM15.all_planets[nr];
                Current.mass *= const_mass/total;
                current_mass += Current.mass;
            }
            cout << "MM15 has mass of " << current_mass << " solar masses." << endl;
        }
        else cout << "MM15 has mass of " << MM15.total_mass << " solar masses." << endl;

        // run system through VV, all data is written to file as we go
        MM15.VelocityVerlet(dimension,integration_points,final_time,force,simple,MM15.total_planets,epsilon);
    }

    return 0;
}



void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << time_step << "; final time = " << final_time << "; integration points = " << N << endl;

    cout << "Initial position = ";
    for(int j=0;j<dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<dimension;j++) cout << v_initial[j] << " ";
    cout << endl;
}

void print_final(int dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}

void randomUniformSphere(double R0,double &x,double &y,double &z,default_random_engine *generator){
    // Random uniform number distribution that returns coordinates (x,y,z) in a sphere of radius R0.

    // Set up the uniform number generator
    uniform_real_distribution<double> uniform_sphere(0.0,1.0);

    // Spherical coordinates
    double phi = 2*M_PI*uniform_sphere(*generator);
    double theta = acos(1 -2*uniform_sphere(*generator));
    double r = R0*pow(uniform_sphere(*generator),1./3);

    // Convert to cartesian coordinates
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);

}

void Gaussian_distribution(double mean,double stddev,double &mass, default_random_engine *generator){
    // Gaussian random number distribution that returns a mass around a mean with a given standard deviation.

    // Set up the uniform number generator
    normal_distribution<double> normal_dist(mean,stddev);

    // Generate the mass
    mass = normal_dist(*generator);
}
