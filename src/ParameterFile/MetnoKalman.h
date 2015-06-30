#ifndef PARAMETER_FILE_METNO_KALMAN_H
#define PARAMETER_FILE_METNO_KALMAN_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"
class ParameterFileMetnoKalman : public ParameterFile {
   public:
      ParameterFileMetnoKalman(const Options& iOptions);

      bool isFixedSize() const {return true;};
      std::vector<int> getTimes() const;

      static bool isValid(std::string iFilename);
      std::string name() const {return "metnoKalman";};
      static std::string description();
   private:
      std::vector<int> mTimes;
      float mLocalMV;
};
#endif
