#ifndef METCAL_SETUP_H
#define METCAL_SETUP_H
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class Calibrator;
class ParameterFile;
class TrainingData;

//! Represents what post-processing method should be trained. Includes which data file to
//! read data from, which output file to place the parameters in, which variable to post-process,
//! and what post-processing method to train.
class SetupTrain {
   public:
      TrainingData* trainingData;
      ParameterFile* output;
      Calibrator* method;
      SetupTrain(const std::vector<std::string>& argv);
      ~SetupTrain();
};
#endif
