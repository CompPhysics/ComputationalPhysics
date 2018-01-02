// normal_distribution example
#include <iostream>
#include <chrono>
#include <random>

void myrandom(long seed){

    std::default_random_engine generator (seed);

  std::normal_distribution<double> distribution (0.0,1.0);

  std::cout << "some Normal-distributed(0.0,1.0) results:" << std::endl;
  for (int i=0; i<10; ++i)
    std::cout << distribution(generator) << std::endl;
}
  int main()
  {
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    myrandom(seed);
    return 0;
  }
