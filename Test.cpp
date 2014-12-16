#include <iostream>
#include <math.h>
#include "Test.h"
#include "assert.h"
#include "DataFile.h"
#include "Variable.h"
#include "Calibration.h"

int main(void) {
   //testDataFile();
   //testWrite();
   testCalibrate();
   ParameterFile pfile("parameters.txt");
   Parameters par = pfile.getParameters(0);
   float ensMean = 0.1519198;
   float ensFrac = 0.9607843;
   int time = 0;
   std::cout << Calibration::getP0(ensMean, ensFrac, par) << std::endl;
   std::cout << Calibration::getP0(0, 1, par) << std::endl;
   for(int i = 1; i <= 9; i++) {
      std::cout << Calibration::getInvCdf(i / 10.0, ensMean, ensFrac, 18, par) << std::endl;
   }
   std::cout << "Testing complete" << std::endl;
   std::cout << mNumErrors << " errors" << std::endl;
}

void testDataFile() {
   std::string filename = "/disk1/eps25_2014010100Z.nc";
   DataFile file(filename);
   std::vector<Variable::Type> types;
   types.push_back(Variable::PrecipAcc);
   types.push_back(Variable::Precip);
   types.push_back(Variable::Cloud);
   for(int v = 0; v < types.size(); v++) {
      const Field& field = file.getField(types[v], 1);
      float total = 0;
      int counter = 0;
      for(int i = 0; i < field.size(); i++) {
         for(int j = 0; j < field[0].size(); j++) {
            for(int e = 0; e < field[0][0].size(); e++) {
               total += field[i][j][e];
               counter++;
            }
         }
      }
      std::cout << "Average " << Variable::getTypeName(types[v]) << ": " << total / counter << std::endl;
   }
}

void testWrite() {
   std::string filename = "test.nc";
   DataFile file(filename);
   Field& field = file.getEmptyField();
   for(int i = 0; i < field.size(); i++) {
      for(int j = 0; j < field[0].size(); j++) {
         for(int e = 0; e < field[0][0].size(); e++) {
            field[i][j][e] = i+j+e;
         }
      }
   }
   file.addField(field,Variable::Precip, 0);
   file.write();
}

void testCalibrate() {
   std::string iFile = "/disk1/eps25_2014010100Z.nc";
   std::string oFile = "test2.nc";
   DataFile fin(iFile);
   DataFile fout(oFile);
   Calibration cal(ParameterFile("test"));
   cal.calibrate(fin, fout);
   fout.write();
}

bool testEq(float iExpected, float iActual, std::string iMessage) {
   float diff = fabs(iExpected - iActual);
   bool isEq = diff / iActual < mFloatTol;
   if(!isEq) {
      std::cout << iMessage << "Value " << iActual << " not equal to the expected value of " << iExpected << std::endl;
      mNumErrors++;
   }
   return isEq;
}
