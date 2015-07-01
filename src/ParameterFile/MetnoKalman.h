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
      //! Applies any translations of values from raw file, such as converting missing values
      //! and potentially applying scale and offset attributes in the future.
      template <class T> void translate(T& value) const {
         if(value == mLocalMV)
            value = Util::MV;
      }
};
#endif
