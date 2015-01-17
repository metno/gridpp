#include "Downscaler.h"

Downscaler::Downscaler() {

}

void Downscaler::downscale(const File& iInput, File& iOutput) const {
   downscaleCore(iInput, iOutput);
}
