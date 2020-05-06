#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../Setup.h"
#include "../KDTree.h"

int main(int argc, const char *argv[]) {
   double start = Util::clock();

   // Create a grid ranging from 30-60 degrees latitude and 30-60 degrees longitude
   int nLat = 1000;
   int nLon = 500;
   std::vector<std::vector<float> > lats;
   std::vector<std::vector<float> > lons;
   lats.resize(nLat);
   lons.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      lats[i].resize(nLon);
      lons[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         lats[i][j] = float(i) / nLat * 30 + 30; // Range 30 to 60
         lons[i][j] = float(j) / nLon * 30 + 30; // Range 30 to 60
      }
   }
   KDTree tree(lats, lons);

   int nRepeat = 100; 

   // Test a point on the inside of the domain
   double s = Util::clock();
   for(int i = 0; i < nRepeat; i++) {
      int I, J;
      tree.getNearestNeighbour(50, 50, I, J);
   }
   double e = Util::clock();
   std::cout << (e - s) / nRepeat << " seconds for points inside domain" << std::endl;

   // Test a point on the outside of the domain
   s = Util::clock();
   for(int i = 0; i < nRepeat; i++) {
      int I, J;
      tree.getNearestNeighbour(10, 10, I, J);
   }
   e = Util::clock();
   std::cout << (e - s) / nRepeat << " seconds for points outside domain" << std::endl;

   return 0;
}
