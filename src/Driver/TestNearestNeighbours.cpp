#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
int main(int argc, const char *argv[]) {
   std::string inputFile     = argv[1];
   std::string outputFile    = argv[2];

   const FileArome ifile(inputFile);
   FileArome ofile(outputFile);
   DownscalerNearestNeighbour downscaler(Variable::T, Options());
   vec2Int I1, I2, J1, J2;
   double t0 = Util::clock();
   downscaler.getNearestNeighbour(ifile, ofile, I1, J1);
   double t1 = Util::clock();
   downscaler.getNearestNeighbourBruteForce(ifile, ofile, I2, J2);
   double t2 = Util::clock();

   std::cout << "Slow: " << t2-t1 << std::endl;
   std::cout << "Fast: " << t1-t2 << std::endl;

   assert(I1.size() == I2.size());
   assert(J1.size() == J2.size());
   for(int i = 0; i < I1.size(); i++) {
      for(int j = 0; j < I1[0].size(); j++) {
         if(I1[i][j] != I2[i][j] || J1[i][j] != J2[i][j]) {
            std::cout << i << " " << j << " " << I1[i][j] << " " << J1[i][j] << " " << I2[i][j] << " " << J2[i][j] << std::endl;
         }
      }
   }

}
