#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;
int main()
{
  vector <vector<int> > vec2D(5, vector<int>(4, 1));
 
  for(auto vec : vec2D)
    {
      for(auto x : vec)
	cout<<x<<" , ";
 
      cout << endl;
    }
 
  cout << endl;
 
  for(int i = 0; i < 5; i++)
    for(int j = 0; j < 5; j++)
      vec2D[i][j] = i*j;
 
  for(auto vec : vec2D)
    {
      for(auto x : vec)
	cout<<x<<" , ";
 
      cout << endl;
    }
 
 
  vec2D.push_back(vector<int>(4, 11));
 
  cout << endl;
 
  for(auto vec : vec2D)
    {
      for(auto x : vec)
	cout<<x<<" , ";
 
      cout << endl;
    }
  return 0;
}






