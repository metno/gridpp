#ifndef SCHEME_H
#define SCHEME_H
class Options;
class ParameterFile;

//! Base class for downscaler, calibrator
class Scheme {
   public:
      Scheme(const Options& iOptions);
      // Is debug turned on? Schemes should display debug information if this is true
      bool debug() const;
   private:
      bool mDebug;
};
#endif
