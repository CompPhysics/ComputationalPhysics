#ifndef VEC3_H
#define VEC3_H

#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

class vec3
{
private:
    double v[3];
//    double *v;

public:
    vec3(); // Default constructor
    ~vec3(); // Destructor

    // Copy constructor
    vec3(vec3 const &copy);

    // Parametetrized constructor
    vec3(double x, double y, double z);

    // Copy assignement constructor
    vec3 &operator= (const vec3 &copy);

    // Printing functions
    void print();

    // Setters
    void setXYZ(double x, double y, double z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    // Getters
    double x() const { return v[0]; }


    // Length

    // Operator overloading
    double &operator[](int index) { return v[index]; }

    // IMPLEMENT MORE OPERATORS HERE
    vec3 &operator+=(const vec3 &other);

    // A cout << foo overload method.
    friend std::ostream &operator<<(std::ostream &os, const vec3 &a) {
        os << "[" << a.v[0] << ", " << a.v[1] << ", " << a.v[2] << "]";
        return os;
    }

};

inline vec3 operator+(vec3 lhs, vec3 rhs)
{
    lhs += rhs;
    return lhs;
}


#endif // VEC3_H
