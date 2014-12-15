#ifndef COEFFS_FILE_H
#define COEFFS_FILE_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "Coeffs.h"

//! Represents a file containing all Fourier coefficients
class CoeffsFile {
   public:
      CoeffsFile(std::string iFilename);
      std::vector<Coeffs> getCoeffsList() const;
   private:
      void parseFile(std::string iFilename);
      std::vector<Coeffs> mCoeffs;
      static bool getNextToken(std::istream& iStream, float& iValue);
};
#endif

