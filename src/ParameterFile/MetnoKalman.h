#ifndef PARAMETER_FILE_METNO_KALMAN_H
#define PARAMETER_FILE_METNO_KALMAN_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"
class ParameterFileMetnoKalman : public ParameterFile {
   public:
      ParameterFileMetnoKalman(std::string iFilename);

      bool isFixedSize() const {return true;};
      std::vector<int> getTimes() const;

      static bool isValid(std::string iFilename);

   private:
      std::vector<int> mTimes;
      float mLocalMV;
};
#endif
