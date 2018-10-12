#include "vec3.h"
#include <cmath>

vec3::vec3()
{
    // Uncomment to see what methods is called when!
//    cout << "Using default constructor" << endl;
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
}


vec3::vec3(double x, double y, double z)
{
    // Uncomment to see what methods is called when!
//    cout << "Using parametrized constructor" << endl;
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

vec3::~vec3()
{
    // Uncomment to see what methods is called when!
//    cout << "Using destructor" << endl;
}

vec3::vec3(vec3 const &copy)
{
    // Uncomment to see what methods is called when!
//    cout << "Using copy constructor" << endl;
    v[0] = copy.v[0];
    v[1] = copy.v[1];
    v[2] = copy.v[2];
}

vec3 &vec3::operator= (const vec3 &copy) {
//    cout << "Using copy assignemnt" << endl;
    v[0] = copy.v[0];
    v[1] = copy.v[1];
    v[2] = copy.v[2];
    return *this;
}

void vec3::print()
{
    cout << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]" << endl;
}


vec3 &vec3::operator+=(const vec3 &other)
{
    v[0] += other.v[0];
    v[1] += other.v[1];
    v[2] += other.v[2];
    return *this;
}
