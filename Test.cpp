#include "Test.h"
#include "assert.h"
#include "Site.h"
#include "Coeffs.h"
#include <iostream>

int main(void) {
   testCoeffs();
   std::cout << "Testing complete" << std::endl;
   std::cout << mNumErrors << " errors" << std::endl;
}

Site getDefaultSite() {
   int id = 1;
   float lat = 61;
   float lon = 5;
   float elev = 230;
   return Site(id, lat, lon, elev);
}

void testCoeffs() {
   Site site = getDefaultSite();
   std::vector<float> coeffsVector;
   coeffsVector.push_back(1);
   coeffsVector.push_back(2);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   coeffsVector.push_back(0);
   mNumErrors = 0;

   const float dirs[]     = {0, 30, 60, 180, 360, 390,  30, 330};
   const float speeds[]   = {1,1,1,1,1,1,              100, 100};
   const float expected[] = {1,2,2.732051,1,1,2,        1.08, 0.92};

   for(int i = 0; i < sizeof(dirs) / sizeof(float); i++) {
      float dir   = dirs[i];
      float speed = speeds[i];
      Coeffs coeffs(site, coeffsVector);
      float factor = coeffs.computeFactor(dir, speed);
      testEq(expected[i], factor);
   }

}

bool testEq(float iExpected, float iActual, std::string iMessage) {
   float diff = fabs(iExpected - iActual);
   bool isEq = diff / iActual < mFloatTol;
   if(!isEq) {
      std::cout << iMessage << "Value " << iActual << " not equal to the expected value of " << iExpected << std::endl;
      mNumErrors++;
   }
   return isEq;
}
