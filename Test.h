#ifndef TEST_H
#define TEST_H
#include <string>

bool testEq(float iExpected, float iActual, std::string iMessage="");
static const float mFloatTol = 1e-5;
static int mNumErrors;
void testDataFile();
void testWrite();
void testCalibrate();
#endif
