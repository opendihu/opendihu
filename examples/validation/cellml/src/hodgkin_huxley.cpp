#include <iostream>
#include <cstdlib>

#include "opendihu.h"

int main(int argc, char *argv[]) {
  // 0D sub-cellular model

  // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv);

  TimeSteppingScheme::Heun<CellmlAdapter<4, 9> // nStates,nAlgebraics: 57,71 =
                                               // Shorten, 4,9 = Hodgkin Huxley
                           >
      equationDiscretized(settings);

  equationDiscretized.run();

  return EXIT_SUCCESS;
}
