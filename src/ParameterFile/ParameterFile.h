#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include <iostream>
#include <map>
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
class ParameterFile {
   public:
      //! Read parameters from this file
      ParameterFile(std::string iFilename);

      //! Get the parameter valid for specified forecast timestep. This is an index, not an hour.
      Parameters getParameters(int iTime) const;

      //! Set the parameter valid for specified time
      void setParameters(Parameters iParameters, int iTime);

      //! Returns the number of parameter sets available (i.e number of forecast hours)
      int getSize() const;

      //! Returns the number of parameters in one parameter set
      int getNumParameters() const;

      //! Returns the filename where parameters are retrieved from
      std::string getFilename() const;
   private:
      std::map<int, Parameters> mParameters; // Offset, Parameters
      std::string mFilename;
      int mNumParameters;
};
#include "ParameterFileSpatial.h"
#endif
