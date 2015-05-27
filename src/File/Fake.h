#ifndef FILE_FAKE_H
#define FILE_FAKE_H
#include "File.h"

//! A file type giving fake data, useful for testing
class FileFake : public File {
   public:
      FileFake(int nLat=10, int nLon=10, int nEns=2, int nTime=10);
      static std::string description();
      std::string name() const {return "fake";};
   protected:
      void writeCore(std::vector<Variable::Type> iVariables);
      FieldPtr getFieldCore(Variable::Type iVariable, int iTime) const;
      bool hasVariableCore(Variable::Type iVariable) const;
};
#endif
