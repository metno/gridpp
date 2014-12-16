#ifndef PARAMETERFILE_H
#define PARAMETERFILE_H
#include <iostream>
#include <map>
#include "Parameters.h"

class ParameterFile {
   public:
      ParameterFile(std::string iFilename);
      Parameters getParameters(int iTime) const;
   private:
      std::map<int, Parameters> mParameters; // Offset, Parameters
      static const int mNumParameters = 9;
      std::string mFilename;
};
#endif
