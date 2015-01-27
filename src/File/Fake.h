#ifndef FILE_FAKE_H
#define FILE_FAKE_H
#include "File.h"

class FileFake : public File {
   public:
      FileFake(int nLat=10, int nLon=10, int nEns=2, int nTime=10);

      int getNumLat() const {return mNLat;};
      int getNumLon() const {return mNLon;};
      int getNumEns() const {return mNEns;};
      int getNumTime() const {return mNTime;};

      vec2 getLats() const;
      vec2 getLons() const;
      vec2 getElevs() const;

      bool setLats(vec2 iLats);
      bool setLons(vec2 iLons);
      bool setElevs(vec2 iElevs);

      std::string name() const {return "fake";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      bool hasVariableCore(Variable::Type iVariable) const {return true;};
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;

      int mNTime;
      int mNLat;
      int mNLon;
      int mNEns;
};
#endif
