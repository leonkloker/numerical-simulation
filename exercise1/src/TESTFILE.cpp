#include "settings.h"

int main(int argc, char *argv[])
{
  Settings test;
  test.loadFromFile(argv[0]);
  test.printSettings();

  return EXIT_SUCCESS;
}
