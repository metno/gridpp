#include "Scheme.h"
#include "Options.h"

Scheme::Scheme(const Options& iOptions) : mDebug(false) {
   iOptions.getValue("debug", mDebug);
}
bool Scheme::debug() const {
   return mDebug;
}
