#include "../Options.h"
#include <gtest/gtest.h>

namespace {
   class OptionsTest : public ::testing::Test {
      protected:
   };

   TEST_F(OptionsTest, accessEmptyOptions) {
      Options options;
      float value = 441;
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(441, value);

      status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(441, value);
   }
   TEST_F(OptionsTest, emptyKey) {
      // Empty keys are not allowed
      Options options;
      std::string value = "new";
      bool status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);

      options.addOption("", "test");
      status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);
   }
   TEST_F(OptionsTest, emptyValue) {
      // Empty values are not allowed
      Options options;
      options.addOption("test", "");
      std::string value = "new";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);
   }
   TEST_F(OptionsTest, emptyValueEqualSign) {
      Options options;
      options.addOptions("test=");
      std::string value = "new";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);

      int valueInt = 2;
      status = options.getValue("test", valueInt);
      EXPECT_FALSE(status);
      EXPECT_EQ(2, valueInt);
   }
   TEST_F(OptionsTest, addOption) {
      Options options;
      options.addOption("test", 441);

      float value = 2;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_FLOAT_EQ(441, value);

      value = 2;
      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(2, value);
   }
   TEST_F(OptionsTest, stringOption) {
      Options options;
      options.addOption("test", "before");

      std::string value = "";
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("before", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("before", value);
   }
   TEST_F(OptionsTest, addOptionEqualSign) {
      Options options;
      options.addOptions("test=before");

      std::string value = "";
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("before", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("before", value);
   }
   TEST_F(OptionsTest, clear) {
      Options options;
      options.addOptions("test=before");
      options.addOptions("other=some");

      options.clear();
      std::string value = "something";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("something", value);

      options.addOptions("new=value");
      status = options.getValue("new", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("value", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("value", value);
   }
   TEST_F(OptionsTest, overwrite) {
      Options options;
      options.addOptions("test=before");
      options.addOptions("test=2");

      float value = 3;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ(2, value);
   }
   TEST_F(OptionsTest, constructor) {
      Options options("test=3 test=2");

      float value = 19;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ(2, value);
      status = options.getValue("new", value);
      EXPECT_FALSE(status);
      EXPECT_EQ(2, value);
   }
   TEST_F(OptionsTest, constructorEmpty) {
      Options options("");

      float value = 19;
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ(19, value);
   }
   TEST_F(OptionsTest, toString) {
      Options options("test=0   q=4");
      std::string s = options.toString();
      EXPECT_LT(0, s.size());

      Options copyOptions(s);
      int i = -1;
      copyOptions.getValue("test", i);
      EXPECT_EQ(0, i);
      copyOptions.getValue("q", i);
      EXPECT_EQ(4, i);

      options.addOptions("f=3");
      s = options.toString();
      copyOptions = Options(s);
      copyOptions.getValue("f", i);
      EXPECT_EQ(3, i);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
