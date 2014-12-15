#ifndef PARAMETERFILE_H
#define PARAMETERFILE_H
#include <iostream>
#include "Parameters.h"

class ParameterFile {
   public:
      ParameterFile(std::string iFilename);
      Parameters getParameters(int iOffset) const;
   private:
      std::vector<Parameters> mParameters; // One set for each offset
};
#endif
