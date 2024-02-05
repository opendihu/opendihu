#pragma once

#include <Python.h> // has to be the first included header
#include <iostream>

namespace Control {

/** A class that measures allocated memory
 */
class MemoryLeakFinder {
public:
  //! get the increase in memory since the last call to this method
  static long long int nKiloBytesIncreaseSinceLastCheck();

  //! output a warning when the memory consumption increases by more than 1 MB
  //! since the last call
  static void warnIfMemoryConsumptionIncreases(std::string message = "");

  //! get the current total memory consumption of the program in kilobytes
  static long long int currentMemoryConsumptionKiloBytes();

private:
  static long long int
      currentMemoryConsumption_; //< the current number of bytes allocated in
                                 //residual set memory
};

} // namespace Control
