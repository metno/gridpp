#include "Cloud.h"
#include "../Util.h"
CalibratorCloud::CalibratorCloud(const ParameterFile& iParameterFile, Variable::Type iPrecip, Variable::Type iCloud) : Calibrator(iParameterFile),
      mCloudType(iCloud),
      mPrecipType(iPrecip) {

}
void CalibratorCloud::calibrateCore(File& iFile) const {
   int nLat = iFile.getNumLat();
   int nLon = iFile.getNumLon();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Initialize calibrated fields in the output file
   std::vector<Field*> cloudCals(nTime);
   for(int t = 0; t < nTime; t++) {
      cloudCals[t] = &iFile.getEmptyField();
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parameters = mParameterFile.getParameters(t);
      const Field& precip = iFile.getField(mPrecipType, t);
      const Field& cloud  = iFile.getField(mCloudType, t);

      Field& cloudCal  = *cloudCals[t];

      // TODO: Figure out which cloudless members to use. Ideally, if more members
      // need precip, we should pick members that already have clouds, so that we minimize
      // our effect on the cloud cover field.

      // Parallelizable
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Turn on clouds if needed, i.e don't allow a member to
            // have precip without cloud cover.
            for(int e = 0; e < nEns; e++) {
               float currPrecip = precip[i][j][e];
               float currCloud  = cloud[i][j][e];
               if(Util::isValid(currPrecip) && Util::isValid(currCloud)) {
                  cloudCal[i][j][e]  = currCloud;
                  // std::cout << "currPrecip = " << currPrecip << std::endl;
                  if(currPrecip > 0 && currCloud < 1) {
                     cloudCal[i][j][e] = 1;
                  }
               }
            }
         }
      }
   }

   // Add to file
   for(int t = 0; t < nTime; t++) {
      iFile.addField(*cloudCals[t], Variable::Cloud, t);
   }
}
