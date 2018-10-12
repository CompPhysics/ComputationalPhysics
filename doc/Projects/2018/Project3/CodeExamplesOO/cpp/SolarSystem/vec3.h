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
    void print(string msg);

    // Setters
    void setXYZ(double x, double y, double z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    // Getters
    double x() const { return v[0]; }
    double y() const { return v[1]; }
    double z() const { return v[2]; }

    // Note: The const after the function name ensures that we
    // cannot modify the returns value, only read it.

    // Length
    double length() const;
    double lengthSquared() const;

    // Operator overloading
    double &operator[](int index) { return v[index]; }


    // IMPLEMENT MORE OPERATORS HERE
    vec3 &operator+=(const vec3 &other);
    vec3 &operator-=(const vec3 &other);

    vec3 &operator*=(const double &other);
    vec3 &operator/=(const double &other);

    // A cout << foo overload method.
    friend std::ostream &operator<<(std::ostream &os, const vec3 &V) {
        os << "[" << V.x() << ", " << V.y() << ", " << V.z() << "]";
        return os;
    }

};

inline vec3 operator+(vec3 lhs, vec3 rhs)
{
    lhs += rhs;
    return lhs;
}

inline vec3 operator-(vec3 lhs, vec3 rhs)
{
    lhs -= rhs;
    return lhs;
}

inline vec3 operator*(vec3 lhs, double rhs)
{
    lhs *= rhs;
    return lhs;
}

inline vec3 operator/(vec3 lhs, double rhs)
{
    lhs /= rhs;
    return lhs;
}


#endif // VEC3_H
