#ifndef TEST_H
#define TEST_H
#include <string>
class Site;

bool testEq(float iExpected, float iActual, std::string iMessage="");
static const float mFloatTol = 1e-5;
static int mNumErrors;
void testCoeffs();
Site getDefaultSite();
#endif
