#ifndef FACTORSFILE_H
#define FACTORSFILE_H
#include "Site.h"
#include <vector>
#include <string>

//! Represents a file where bias-factors are stored
class FactorsFile {
   public:
      FactorsFile(int iDate, int iInit);
      //! Add the factors to the queue
      void add(Site iSite, float iModelElev, std::vector<float> iFactors);
      //! Write all queued facotors to specified filename
      void write(std::string iFilename);
   private:
      std::vector<Site> mSites;
      std::vector<std::vector<float> > mFactors;
      std::vector<float> mModelElevs;
      int mDate;
      int mInit;
};
#endif

