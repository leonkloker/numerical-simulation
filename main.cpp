#include <iostream>
#include <cstdlib>
#include "settings.h"

int main(int argc, char *argv[])
{
  #ifndef NDEBUG
  // only run this code in debug target
  std::cout << "lots of inefficient but informative output . . .";
  #endif

  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];
  
  // print message
  std::cout << "Filename: \"" << filename << "\"" << std::endl;

  loadFromFile(filename);

  return EXIT_SUCCESS;
}
