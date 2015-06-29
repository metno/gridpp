#ifndef PARAMETER_FILE_TEXT_H
#define PARAMETER_FILE_TEXT_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"

//! Represents a collection of parameters, one set for each forecast time
//! Parameters are read from a specified text file with the following format:
//! 0 3.4 2.1 5.2 12 41
//! 1 3.4 2.1 5.2 12 41
//! 2 3.4 2.1 5.2 12 41
//! Each line represents one forecast time. The first column is the forecast timestep (an index
//! not the number of hours), starting at 0. The remaning columns are parameters that can be used
//! in post-processing methods. The number of columns most be constant. If a file has only one
//! line, then the parameters are used for all forecast hours.
class ParameterFileText : public ParameterFile {
   public:
      //! Read parameters from this file
      ParameterFileText(std::string iFilename, bool iIsSpatial=false);

      bool isFixedSize() const;

      std::vector<int> getTimes() const;
      static bool isValid(std::string iFilename);

      int getNumParameters() const;
      //! Write parameter file to disk
      //! @param iFilename Write to this filename. If empty, write to the file that was read from.
      void write(const std::string& iFilename="") const;
   private:
      int mNumParameters;
      bool mIsSpatial;
      std::vector<int> mTimes;
};
#endif
