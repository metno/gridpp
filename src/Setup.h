#ifndef METCAL_SETUP_H
#define METCAL_SETUP_H
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class File;
class Calibrator;
class Downscaler;

//! Represents the post-processing of one variable
struct VariableConfiguration {
   //! Which variable should be post-processed?
   Variable::Type variable;
   //! Which downscaler should be used
   Downscaler* downscaler;
   //! Which calibrators should be use
   std::vector<Calibrator*> calibrators;
   Options variableOptions;
};

//! Represents what and how the post-processing should be done. Includes which input file to
//! post-process, which output file to place the results in, which variables to post-process,
//! and what post-processing methods to invoke on each variable.
class Setup {
   public:
      File* inputFile;
      File* outputFile;
      std::vector<VariableConfiguration> variableConfigurations;
      Setup(std::vector<std::string> argv);
      ~Setup();
      static std::string defaultDownscaler();
};
#endif
