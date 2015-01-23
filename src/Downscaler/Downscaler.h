#ifndef DOWNSCALER_H
#define DOWNSCALER_H
#include <string>
#include "../Options.h"
#include "../Variable.h"
class File;

class Downscaler {
   public:
      Downscaler(Variable::Type iVariable);
      void downscale(const File& iInput, File& iOutput) const;
      static Downscaler* getScheme(std::string iName, Variable::Type iVariable, Options& iOptions);
      virtual std::string name() const = 0;
   protected:
      virtual void downscaleCore(const File& iInput, File& iOutput) const = 0;
      Variable::Type mVariable;
};
#include "NearestNeighbour.h"
#include "Gradient.h"
#include "Smart.h"
#endif
