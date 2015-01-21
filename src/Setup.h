#ifndef METCAL_SETUP_H
#define METCAL_SETUP_H
#include <string>
#include <vector>
#include "Variable.h"
class File;
class Calibrator;
class Downscaler;

struct VariableConfiguration {
   Variable::Type variable;
   Downscaler* downscaler;
   std::vector<Calibrator*> calibrators;
};
class Setup {
   public:
      File* inputFile;
      File* outputFile;
      std::vector<VariableConfiguration> variableConfigurations;
      static bool getSetup(std::vector<std::string> argv, Setup& iSetup);
};
#endif
