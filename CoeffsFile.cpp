#include "CoeffsFile.h"
#include <assert.h>
#include <stdlib.h>

CoeffsFile::CoeffsFile(std::string iFilename) {
   parseFile(iFilename);
}

std::vector<Coeffs> CoeffsFile::getCoeffsList() const {
   return mCoeffs;
}

void CoeffsFile::parseFile(std::string iFilename) {

   std::ifstream ifs(iFilename.c_str(), std::ifstream::in);
   assert(ifs.good());
   char line[10000]; 
   while(ifs.good()) {
      ifs.getline(line, 10000, '\n');   
      if(ifs.good() && line[0] != '#') {   
         std::stringstream ss(line);
         float id;
         float lat;
         float lon;
         float elev;
         getNextToken(ss, id);
         getNextToken(ss, lat);
         getNextToken(ss, lon);
         getNextToken(ss, elev);
         std::vector<float> coeffs;
         while(true) {
            float coeff;
            bool status = getNextToken(ss, coeff);
            if(status == false)
               break;
            coeffs.push_back(coeff);

         }
         Site site((int) id, lat, lon, elev);
         mCoeffs.push_back(Coeffs(site, coeffs));
      }
   }
}
bool CoeffsFile::getNextToken(std::istream& iStream, float& iValue) {
   std::string str;
   std::getline(iStream, str, ',');
   if(str == "")
      return false;
   if(str == "NA") {
      std::cout << "ERROR: Coefficients file has NAs" << std::endl;
      abort();
   }
   std::stringstream ss(str);
   ss >> iValue;
   return true;
}
