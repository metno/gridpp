#ifndef METCAL_SETUP_H
#define METCAL_SETUP_H
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class File;
class Calibrator;
class Downscaler;

struct VariableConfiguration {
   Variable::Type variable;
   Downscaler* downscaler;
   std::vector<Calibrator*> calibrators;
   Options variableOptions;
};
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
