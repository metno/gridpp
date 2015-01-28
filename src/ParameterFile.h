#ifndef PARAMETERFILE_H
#define PARAMETERFILE_H
#include <iostream>
#include <map>
#include "Parameters.h"

class ParameterFile {
   public:
      ParameterFile(std::string iFilename);
      Parameters getParameters(int iTime) const;
      void       setParameters(Parameters iParameters, int iTime);
      int        getSize() const;
      int        getNumParameters() const;
      std::string getFilename() const;
   private:
      std::map<int, Parameters> mParameters; // Offset, Parameters
      std::string mFilename;
      int mNumParameters;
};
/*
class ParameterFileRegion {
   public:
      ParameterFileRegion(std::string iFilename) {};
      Parameters getParameters(float iLat, float iLon, int iTime) const {return Parameters();};
   private:
      std::map<int, Parameters> mParameters; // Offset, Parameters
      std::string mFilename;
};
*/
#endif
