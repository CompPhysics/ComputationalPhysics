#include <iostream>
#include "vec3.h"

using namespace std;

int main()
{
    cout << "   This shows how to use our new 3-vector class. We can create a vector like this: " << endl;
    cout << "vec3 myVector(1,2,3);" << endl;
    vec3 myVector(1,2,3);
    cout << "   We can print it like this: " << endl;
    cout << "myVector.print();" << endl;
    myVector.print();
    cout << "   Or even with a name for fast copy/paste into Matlab" << endl;
    cout << "myVector.print(\"v\");" << endl;
    myVector.print("v");
    cout << "   Another way to print the vector" << endl;
    cout << "cout << myVector << endl;" << endl;
    cout << myVector << endl;
    cout << "   You can get the length (or norm) of the vector like this" << endl;
    cout << "myVector.length()" << endl;
    cout << myVector.length() << endl;

    cout << "   We can also access elements with either square brackets or normal parentheses" << endl;
    cout << "myVector[0]" << endl;
    cout << myVector[0] << endl;
    cout << "myVector(1)" << endl;
    cout << myVector(1) << endl;
    cout << "   Or use" << endl;
    cout << "myVector.x()" << endl;
    cout << myVector.x() << endl;
    cout << "   You can set elements with the first two as well" << endl;
    myVector[0] = 5;
    myVector.print("v");

    cout << "   We can now also do arithmetic on the vector. If you want to add 1 to all components, we can write" << endl;
    cout << "myVector += 1;" << endl;
    myVector += 1;
    myVector.print("v");
    cout << "   Or multiply by 2.3" << endl;
    cout << "myVector *= 10.0;" << endl;
    myVector *= 10.0;
    myVector.print("v");
    cout << "   We can even do arithmetic with a left hand side scalar " << endl;
    cout << "2+myVector" << endl;
    cout << 2+myVector << endl;

    cout << "   What about arithmetic on two vectors. Let's create two vectors a and b. But this time with a new syntax (notice the curly brackets)" << endl;
    cout << "vec3 a({1,2,3});" << endl;
    cout << "vec3 b({2,3,4});" << endl;
    vec3 a({1,2,3});
    vec3 b({2,3,4});
    cout << "   The sum of the two is" << endl;
    cout << "cout << a+b << endl;" << endl;
    cout << a+b << endl;
    cout << "   which has the length" << endl;
    cout << "cout << (a+b).length() << endl;" << endl;
    cout << (a+b).length() << endl;

    cout << "   The cross product between the two is" << endl;
    cout << "a.cross(b); " << endl;
    cout << a.cross(b) << endl;

    cout << "   We can also create a new vector and set it to an existing one" << endl;
    cout << "vec3 c = a;" << endl;
    vec3 c = a;
    c.print("c");

    return 0;
}










