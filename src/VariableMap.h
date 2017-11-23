#ifndef VARIABLE_MAP_H
#define VARIABLE_MAP_H
#include <map>
#include <string>
#include "Variable.h"

class VariableMap {
   public:
      VariableMap();
      void add(Variable::Type iType, std::string iName);
      bool has(Variable::Type iType) const;
      //! Get the name of a variable type
      std::string get(Variable::Type iType) const;
      static std::vector<Variable> load(std::string iFilename);
      static std::vector<Variable> loadDefaults();
      static Variable getDefault(Variable::Type);
      void clear();
   private:
      std::map<int, std::string> mMap; // Variable::Type, variable name in file
};
#endif

