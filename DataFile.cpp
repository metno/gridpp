#include "DataFile.h"
#include <math.h>
#include <assert.h>

DataFile::DataFile(std::string iFilename) :
   mFilename(iFilename), mFile(iFilename.c_str()) {
   if(!mFile.is_valid()) {
      std::cout << "Error: Netcdf file " << iFilename << " not valid" << std::endl;
   }
}

void DataFile::read() {
   // Set offsets
   NcDim* dimTime = mFile.get_dim("time");
   mNumOffsets = dimTime->size();
   // Set lat/lons
   int numX = mFile.get_dim("x")->size(); 
   int numY = mFile.get_dim("y")->size(); 
   NcVar* varLats = mFile.get_var("latitude"); 
   NcVar* varLons = mFile.get_var("longitude"); 
   NcVar* varElevs = mFile.get_var("surface_geopotential"); 

   float* lats = new float[numX*numY];
   float* lons = new float[numX*numY];
   float* elevs = new float[numX*numY];
   long count2D[2] = {numY, numX};
   long count4D[4] = {1, 1, numY, numX};
   varLats->get(lats, count2D);
   varLons->get(lons, count2D);
   varElevs->get(elevs, count4D);

   mLats.resize(numX*numY);
   mLons.resize(numX*numY);
   mElevs.resize(numX*numY);

   int index = 0;
   for(int x = 0; x < numX; x++) {
      for(int y = 0; y < numY; y++) {
         mLats[index] = lats[index];
         mLons[index] = lons[index];
         mElevs[index] = elevs[index] / 9.81;
         index++;
      }
   }

   delete[] lats;
   delete[] lons;
   delete[] elevs;
}

const Field& DataFile::getField(std::string iVariable, int iTime) const {
   std::map<std::string, std::vector<Field> >::const_iterator it = mValues.find(iVariable);
   if(it == mValues.end()) {
      // Not cached, retrieve data
      NcVar* var = mFile.get_var(iVariable.c_str());
      int numT = getNumOffsets();
      int numE = 1;
      int numLat = 1;
      int numLon = 1;
      long count[5] = {numT, 1, numE, numLat, numLon};
      float* values = new float[numT*1*numE*numLat*numLon];
      var->get(values, count);
      int index = 0;
      for(int t = 0; t < numT; t++) {
         Field field;
         for(int lat = 0; lat < numLat; lat++) {
            for(int lon = 0; lon < numLon; lon++) {
               for(int e = 0; e < numE; e++) {
                  field[lat][lon][e] = values[index];
                  index++;
               }
            }
         }
         mValues[iVariable].push_back(field);
      }
      delete[] values;
   }
   return it->second[iTime];
}

DataFile::~DataFile() {
   mFile.close();
}


int DataFile::getDate() const {
   return mDate;
}
int DataFile::getInit() const {
   return mInit;
}
