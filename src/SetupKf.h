#ifndef METCAL_SETUP_KF_H
#define METCAL_SETUP_KF_H
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class Calibrator;
class Downscaler;
class ParameterFile;
class File;

//! Represents how the Kalman filter should be run. Includes which data file to
//! read data from, which kalman database files to read and update, and which output file to place
//! the estimated biases in.
class SetupKf {
   public:
      ParameterFile* output;
      File* fcstFile;
      File* obsFile;
      Calibrator* method;
      Downscaler* downscaler; // Not used
      Variable::Type variable;
      ParameterFile* dbin;
      ParameterFile* dbout;
      Options options;
      SetupKf(const std::vector<std::string>& argv);
      ~SetupKf();
};
#endif
