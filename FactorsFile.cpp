#include "FactorsFile.h"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>

FactorsFile::FactorsFile(int iDate, int iInit) : mDate(iDate), mInit(iInit) {

}
void FactorsFile::add(Site iSite, float iModelElev, std::vector<float> iFactors) {
   mSites.push_back(iSite);
   mModelElevs.push_back(iModelElev);
   mFactors.push_back(iFactors);
}
void FactorsFile::write(std::string iFilename) {
   assert(mSites.size() == mFactors.size());
   std::ofstream ofs(iFilename.c_str());

   int year = (mDate / 10000);
   int month = (mDate / 100) % 100;
   int day = mDate % 100;
   int hour = mInit;
   int nobs = mSites.size();
   if(nobs == 0) {
      std::cout << "Warning: No observations available. No file written" << std::endl;
      return;
   }
   int nhours = mFactors[0].size();

   // Write header
   ofs << year << " " << month << " " << day << " " << hour << std::endl;
   ofs << nobs << " " << nhours << std::endl;
   ofs << std::fixed;
   for(int i = 0; i < mSites.size(); i++) {
      float modelHeight = mModelElevs[i];
      std::vector<float> factors = mFactors[i];
      assert(factors.size() == nhours);
      Site site = mSites[i];
      ofs << std::setprecision(0)
          << std::setw(6) << site.getId()
          << std::setprecision(4)
          << std::setw(9) << site.getLat()
          << std::setw(9) << site.getLon()
          << std::setprecision(0)
          << std::setw(9) << site.getElev()
          << std::setprecision(2)
          << std::setw(9) << modelHeight;
      for(int j = 0; j < factors.size(); j++) {
         ofs << std::setw(7) << std::setprecision(2) << factors[j];
      }
      ofs << std::endl;
   }
   ofs.close();
}
