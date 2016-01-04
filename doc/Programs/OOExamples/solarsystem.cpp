#include "solarsystem.h"
#include "planet.h"
#include <iostream>
#include <armadillo>
#include <stdio.h>
#include <iomanip>
using namespace arma;
using namespace std;


solarsystem::solarsystem()
{
}
void solarsystem::add(planet n){
     number_planets++;
     all_planets.push_back(n);
 }
void solarsystem::print_position(ofstream &output, ofstream &output2, vector<planet> vec){
    print_position(output, output2, vec, 3);
}

void solarsystem::print_position(ofstream &output, ofstream &output2, vector<planet> vec, int n){
    if(n>3 || n<=0) n=3;
    for(int i=0; i<vec.size(); i++){
        planet &this = vec[i];
        std::cout << std::scientific;
        for(int j=0; j<n;j++){
        std::cout << this.position[j] << "   ";
        output << std::scientific << this.position[j] << "   ";
        output2 << std::scientific << this.velocity[j] << "   ";
        }
        std::cout << "         ";
        output  << "         ";
        output2  << "         ";
       }
    std::cout << std::endl;
    output << endl;
    output2 << endl;
}

void solarsystem::insert_data(vector<planet> vec, mat &ma){
    for(int i=0; i<vec.size(); i++){
        planet &this = vec[i];
        ma(i,6)=this.mass;

        for(int k=0; k<3;k++){
            ma(i,k)=this.position[k];
            ma(i,k+3)=this.velocity[k];
        }
    }
}

void solarsystem::synctroniz(vector<planet> vec, mat &ma){
    int n = vec.size();

    for(int j=0; j<n;j++){
        planet &this = vec[j];
    for (int i = 0; i < 3; ++i){
       this.position[i] =  ma(j,i);
       this.velocity[i] = ma(j,i+3);
}
}

}

void solarsystem::solverRK4(vector<planet> vec, double h, double tmax){

    mat y_i(number_planets,7);
    mat y_i_temp(number_planets,7);
    mat k1(number_planets,7);
    mat k2(number_planets,7);
    mat k3(number_planets,7);
    mat k4(number_planets,7);

    insert_data(vec , y_i);
    double t=0;

    //for print file_________
    char *filename = new char[1000];
    char *filename2 = new char[1000];
        sprintf(filename, "Planet_position_RK4_%f.dat", h);
        sprintf(filename2, "Planet_velocity_RK4_%f.dat", h);

        ofstream output (filename);
        ofstream output2 (filename2);

        if (output.is_open()){
            output.precision(5);
            output2.precision(5);
    while(t<tmax){

        derivative(y_i, k1, number_planets);

        sum_matrix(y_i_temp, 1, y_i, 0.5*h, k1, number_planets);
        derivative(y_i_temp, k2, number_planets);

        sum_matrix( y_i_temp, 1,  y_i, 0.5*h,  k2, number_planets);

        derivative( y_i_temp,  k3, number_planets);

        sum_matrix( y_i_temp, 1,  y_i, h,  k3, number_planets);

        derivative( y_i_temp,  k4, number_planets);

        for(int j=1; j<number_planets; j++){

             for(int i=0; i<6; i++){
                 y_i(j,i) = y_i(j,i) + h*(k1(j,i) + 2*k2(j,i) + 2*k3(j,i) + k4(j,i))/6;
             }

             //Syncroniz position and velocity with the classes
             planet &this = vec[j];
             for(int i=0; i<3; i++){
             this.position[i] = y_i(j,i);
             this.velocity[i] = y_i(j,i+3);
             }

        }

print_position(output, output2, vec, 3);

t+=h;

}

output.close();
}
}

void solarsystem::solverVERLET(vector<planet> vec, double h, double tmax){

    mat y_i(number_planets,7);
    mat r_i_dt(number_planets,7);
    mat a_dt(number_planets,7);
    mat v_dt(number_planets,7);
    mat v_dt_2(number_planets,7);

    insert_data(vec , y_i);

    double t=0,zz =1;

    //print file___________________
    char *filename = new char[1000];
    char *filename2 = new char[1000];
        sprintf(filename, "Planet_position_Verlet_%f.dat", h);
        sprintf(filename2, "Planet_velocity_Verlet_%f.dat", h);

            ofstream output (filename);
            ofstream output2 (filename2);

            if (output.is_open()){
            output.precision(5);
            output2.precision(5);
    // end for print


    while(t<tmax){

        derivative(y_i, a_dt, number_planets);

        for(int j=0; j<number_planets; j++){

             for(int i=0; i<3; i++){
                 y_i(j,i) = y_i(j,i) + h*y_i(j,i+3) + 0.5*h*h*a_dt(j,i+3);
                 v_dt_2(j,i+3) = y_i(j,i+3) + 0.5*h*a_dt(j,i+3);
             }
        }

        derivative(y_i, a_dt, number_planets);

        for(int j=0; j<number_planets; j++){

             for(int i=3; i<6; i++){

                 y_i(j,i) = v_dt_2(j,i) + 0.5*h*a_dt(j,i);

             }

             planet &this = vec[j];
             for(int i=0; i<3; i++){
             this.position[i] = y_i(j,i);
             this.velocity[i] = y_i(j,i+3);
             }
        }
print_position(output, output2, vec,3);

t+=h;

}
output.close();
}


}
void solarsystem::sum_matrix(mat &result, double coeff_one, mat &first,double coeff_two, mat &second, int n){
    for(int j=0; j<n; j++){
         for(int i=0; i<6; i++){
            result(j,i) = coeff_one*first(j,i) + coeff_two*second(j,i);
         }
         result(j,6) = first(j,6);

    }
}
void solarsystem::printmat(mat &ma, int n){
    cout << endl;
    for(int i=0; i<7; i++){

        for(int k=0; k<n;k++){
            cout <<  ma(k,i)<<" " ;
        } cout << endl;}
}
double solarsystem::force(double x, double y, double z, double Mothers){
    double G=  4*M_PI*M_PI;
    double force=0;
    double distance=0;

    distance = x*x + y*y + z*z;

    force = G*Mothers/pow(distance, 1.5);

    return force;
}
void solarsystem::derivative(mat &dat, mat &de, int n){

    double accelleration_x=0,accelleration_y=0,accelleration_z=0, mod_force;
    for(int i=0; i<n; i++){

        accelleration_x=0,accelleration_y=0,accelleration_z=0;
        for(int j=0; j<n; j++){
            if(i!=j){

mod_force = force(dat(j,0)-dat(i,0),dat(j,1)-dat(i,1) ,dat(j,2)-dat(i,2),  dat(j,6));


accelleration_x += mod_force*(dat(j,0)-dat(i,0));
accelleration_y += mod_force*(dat(j,1)-dat(i,1));
accelleration_z += mod_force*(dat(j,2)-dat(i,2));
}
}
        de(i,3) = accelleration_x;
        de(i,4) = accelleration_y;
        de(i,5) = accelleration_z;

    }


    for(int i=0; i<n; i++){
        de(i,0) = dat(i,3); //velx
        de(i,1) = dat(i,4); //vely
        de(i,2) = dat(i,5); //velz
    }
}
