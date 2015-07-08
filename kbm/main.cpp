#include <iostream>

#include <gsl/gsl_combination.h>

#include "kSystem.h"
#include "kCluster.h"
#include "kMath.h"

using namespace std;

int main(int argc, const char *argv[]) {
  std::cout << std::endl << "Started...\n";
  kSystem *mysystem = new kSystem();

  delete (mysystem);

  std::cout << std::endl << "\n...Finished!\n\n";
  return 0;
}
