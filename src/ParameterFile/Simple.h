#ifndef PARAMETER_FILE_SIMPLE_H
#define PARAMETER_FILE_SIMPLE_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"

class ParameterFileSimple : public ParameterFile {
   public:
      //! Read parameters from this file
      ParameterFileSimple(Parameters iParameters);

      bool isFixedSize() const {return true;};

      std::vector<int> getTimes() const;
      static bool isValid(std::string iFilename);
      std::string name() const {return "simple";};
   private:
      int mNumParameters;
      std::vector<int> mTimes;
};
#endif
