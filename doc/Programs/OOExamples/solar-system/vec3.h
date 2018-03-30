#ifndef VEC3_H
#define VEC3_H
#include <string>
#include <vector>

class vec3
{
public:
    vec3();
    vec3(double x, double y, double z);
    double lengthSquared();
    double length();

    // Actions
    void zeros();
    vec3 cross(vec3 otherVector);
    double dot(vec3 otherVector);
    void normalize();
    vec3 normalized();

    // Getters and setters
    double x() const { return components[0]; }
    double y() const { return components[1]; }
    double z() const { return components[2]; }
    void setX(double x) { components[0] = x; }
    void setY(double y) { components[1] = y; }
    void setZ(double z) { components[2] = z; }

    // Convenience functions
    void print();
    void print(std::string name);
    friend std::ostream& operator<<(std::ostream& os, const vec3& myVector); // Allows cout << myVector << endl;

    // Operators
    double &operator()(int index) { return components[index]; } // Allows access like myVector(0)
    double &operator[](int index) { return components[index]; } // Allows access like myVector[0]
    vec3 &operator+=(double rhs); // Componentwise addition with scalar
    vec3 &operator+=(vec3 rhs);   // Componentwise addition with vector
    vec3 &operator*=(double rhs); // Componentwise multiplication with scalar
    vec3 &operator*=(vec3 rhs);   // Componentwise multiplicationwith vector
    vec3 &operator-=(double rhs); // Componentwise subtraction with scalar
    vec3 &operator-=(vec3 rhs);   // Componentwise subtraction with vector
    vec3 &operator/=(double rhs); // Componentwise division with scalar
    vec3 &operator/=(vec3 rhs);   // Componentwise division with vector
private:
    double components[3];
};

inline vec3 operator+(vec3 lhs, double rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator+(double lhs, vec3 rhs) {
    rhs += lhs;
    return rhs;
}

inline vec3 operator+(vec3 lhs, vec3 rhs) {
    lhs += rhs;
    return lhs;
}

inline vec3 operator-(vec3 lhs, double rhs) {
    lhs -= rhs;
    return lhs;
}

inline vec3 operator-(double lhs, vec3 rhs) {
    rhs -= lhs;
    return rhs;
}

inline vec3 operator-(vec3 lhs, vec3 rhs) {
    lhs -= rhs;
    return lhs;
}


inline vec3 operator*(vec3 lhs, double rhs) {
    lhs *= rhs;
    return lhs;
}

inline vec3 operator*(double lhs, vec3 rhs) {
    rhs *= lhs;
    return rhs;
}

inline vec3 operator*(vec3 lhs, vec3 rhs) {
    lhs *= rhs;
    return lhs;
}


inline vec3 operator/(vec3 lhs, double rhs) {
    lhs /= rhs;
    return lhs;
}

inline vec3 operator/(double lhs, vec3 rhs) {
    rhs /= lhs;
    return rhs;
}

inline vec3 operator/(vec3 lhs, vec3 rhs) {
    lhs /= rhs;
    return lhs;
}

#endif // VEC3_H
