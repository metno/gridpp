#include "gridpp.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <assert.h>
#include <execinfo.h>
#include <signal.h>
#include <iomanip>
#include <cstdio>

#ifdef DEBUG
extern "C" void __gcov_flush();
#endif

bool gridpp::util::is_valid(float value) {
    return !std::isnan(value) && !std::isinf(value) && value != gridpp::MV;
}
float gridpp::util::calculate_stat(const std::vector<float>& iArray, gridpp::util::StatType iStatType, float iQuantile) {
   // Initialize to missing
   float value = gridpp::MV;
   if(iStatType == gridpp::util::StatTypeMean || iStatType == gridpp::util::StatTypeSum) {
      float total = 0;
      int count = 0;
      for(int n = 0; n < iArray.size(); n++) {
         if(gridpp::util::is_valid(iArray[n])) {
            total += iArray[n];
            count++;
         }
      }
      if(count > 0) {
         if(iStatType == gridpp::util::StatTypeMean)
            value = total / count;
         else
            value = total;
      }
   }
   else if(iStatType == gridpp::util::StatTypeStd) {
      // STD = sqrt(E[X^2] - E[X]^2)
      // The above formula is unstable when the variance is small and the mean is large.
      // Use the property that VAR(X) = VAR(X-K). Provided K is any element in the array,
      // the resulting calculation of VAR(X-K) is stable. Set K to the first non-missing value.
      float total  = 0;
      float total2 = 0;
      float K = gridpp::MV;
      int count = 0;
      for(int n = 0; n < iArray.size(); n++) {
         if(gridpp::util::is_valid(iArray[n])) {
            if(!gridpp::util::is_valid(K))
               K = iArray[n];
            assert(gridpp::util::is_valid(K));

            total  += iArray[n] - K;
            total2 += (iArray[n] - K)*(iArray[n] - K);
            count++;
         }
      }
      if(count > 0) {
         float mean  = total / count;
         float mean2 = total2 / count;
         float var   = mean2 - mean*mean;
         if(var < 0) {
            // This should never happen
            var = 0;
            // Util::warning("CalibratorNeighbourhood: Problems computing std, unstable result. Setting value to 0");
         }
         float std = sqrt(var);
         value = std;
      }
   }
   else {
      if(iStatType == gridpp::util::StatTypeMin)
         iQuantile = 0;
      if(iStatType == gridpp::util::StatTypeMedian)
         iQuantile = 0.5;
      if(iStatType == gridpp::util::StatTypeMax)
         iQuantile = 1;
      // Remove missing
      std::vector<float> cleanHood;
      cleanHood.reserve(iArray.size());
      for(int i = 0; i < iArray.size(); i++) {
         if(gridpp::util::is_valid(iArray[i]))
            cleanHood.push_back(iArray[i]);
      }
      int N = cleanHood.size();
      if(N > 0) {
         std::sort(cleanHood.begin(), cleanHood.end());
         int lowerIndex = floor(iQuantile * (N-1));
         int upperIndex = ceil(iQuantile * (N-1));
         float lowerQuantile = (float) lowerIndex / (N-1);
         float upperQuantile = (float) upperIndex / (N-1);
         float lowerValue = cleanHood[lowerIndex];
         float upperValue = cleanHood[upperIndex];
         if(lowerIndex == upperIndex) {
            value = lowerValue;
         }
         else {
            assert(upperQuantile > lowerQuantile);
            assert(iQuantile >= lowerQuantile);
            float f = (iQuantile - lowerQuantile)/(upperQuantile - lowerQuantile);
            assert(f >= 0);
            assert(f <= 1);
            value   = lowerValue + (upperValue - lowerValue) * f;
         }
      }
   }
   return value;
}
int gridpp::util::num_missing_values(const vec2& iArray) {
   int count = 0;
   for(int y = 0; y < iArray.size(); y++) {
      for(int x = 0; x < iArray[y].size(); x++) {
         count += !gridpp::util::is_valid(iArray[y][x]);
      }
   }
   return count;
}
void gridpp::util::debug(std::string string) {
    std::cout << string << std::endl;
}

void gridpp::util::error(std::string iMessage) {
#ifdef DEBUG
  std::cout << "Error: " << iMessage << std::endl;
  void *array[10];
  size_t size = backtrace(array, 10);
  std::cout << "Stack trace:" << std::endl;
  backtrace_symbols_fd(array, size, 2);
   __gcov_flush();
#else
  std::cout << "Error: " << iMessage << std::endl;
#endif
   abort();
}
double gridpp::util::clock() {
   timeval t;
   gettimeofday(&t, NULL);
   double sec = (t.tv_sec);
   double msec= (t.tv_usec);
   return sec + msec/1e6;
}
gridpp::util::StatType gridpp::util::getStatType(std::string iName) {
    gridpp::util::StatType iType;
   if(iName == "mean") {
      iType = gridpp::util::StatTypeMean;
   }
   else if(iName == "min") {
      iType = gridpp::util::StatTypeMin;
   }
   else if(iName == "max") {
      iType = gridpp::util::StatTypeMax;
   }
   else if(iName == "median") {
      iType = gridpp::util::StatTypeMedian;
   }
   else if(iName == "quantile") {
      iType = gridpp::util::StatTypeQuantile;
   }
   else if(iName == "std") {
      iType = gridpp::util::StatTypeStd;
   }
   else if(iName == "sum") {
      iType = gridpp::util::StatTypeSum;
   }
   return iType;
}
vec gridpp::util::calc_even_quantiles(const vec& values, int num) {
    if(num >= values.size())
        return  values;

    vec quantiles;
    if(num == 0)
        return quantiles;

    quantiles.reserve(num);
    quantiles.push_back(values[0]);
    int count_lower = 0;
    for(int i = 0; i < values.size(); i++) {
        if(values[i] != quantiles[0])
            break;
        count_lower++;
    }
    int size = values.size();
    if(count_lower > size / num)
        quantiles.push_back(values[count_lower]);

    // Remove duplicates
    vec values_unique;
    values_unique.reserve(values.size());
    float last_quantile = quantiles[quantiles.size() - 1];
    for(int i = 0; i < values.size(); i++) {
        if(values[i] > last_quantile && (values_unique.size() == 0 || values[i] != values_unique[values_unique.size() - 1]))
            values_unique.push_back(values[i]);
    }
    int num_left = num - quantiles.size();
    // std::cout << "Number of unique values: " << values_unique.size() << std::endl;
    for(int i = 1; i <= num_left; i++) {
        float f = float(i) / (num_left);
        int index = values_unique.size() * f - 1;
        if(index > 0) {
            float value = values_unique[index];
            quantiles.push_back(value);
        }
        else {
            std::cout << i << " " << f << " " << index << " " << num_left << " " << values_unique.size() << std::endl;
            std::cout << count_lower << " " << values.size() << " " << last_quantile << std::endl;
            abort();
        }
    }
    return quantiles;
    // for(int i = 0; i < quantiles.size(); i++)
    //     std::cout << "Threshold[" << i << "] = " << quantiles[i] << std::endl;
    /*
    int index = size / num;
    int step = size / num;
    count = 1;
    while(index < Y * X * E) {
        float value = values[index];
        std::cout << " " << count << " " << index << " " << step << " " << quantiles.size() << std::endl;
        if(value != last) {
            std::cout << "Threshold " << quantiles.size() << " = " << value << std::endl;
            quantiles.push_back(value);
            last = value;
        }
        else {
            step = (size - index) / (num - quantiles.size());
            std::cout << "Same. Step is now: " << step << std::endl;
        }
        last = value;
        index += step;
        count++;
    }
    */

}
int gridpp::util::get_lower_index(float iX, const std::vector<float>& iValues) {
    int index = gridpp::MV;
    for(int i = 0; i < (int) iValues.size(); i++) {
        float currValue = iValues[i];
        if(gridpp::util::is_valid(currValue)) {
            if(currValue < iX) {
                index = i;
            }
            else if(currValue == iX) {
                index = i;
                break;
            }
            else if(currValue > iX) {
                break;
            }
        }
    }
    return index;
}
int gridpp::util::get_upper_index(float iX, const std::vector<float>& iValues) {
    int index = gridpp::MV;
    for(int i = iValues.size()-1; i >= 0; i--) {
        float currValue = iValues[i];
        if(gridpp::util::is_valid(currValue)) {
            if(currValue > iX) {
                index = i;
            }
            else if(currValue == iX) {
                index = i;
                break;
            }
            else if(currValue < iX) {
                break;
            }
        }
    }
    return index;
}
float gridpp::util::interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY) {
    float y = gridpp::MV;

    if(x >= iX[iX.size()-1])
        return iY[iX.size()-1];
    if(x <= iX[0])
        return iY[0];

    int i0   = get_lower_index(x, iX);
    int i1   = get_upper_index(x, iX);
    float x0 = iX[i0];
    float x1 = iX[i1];
    float y0 = iY[i0];
    float y1 = iY[i1];

    if(x0 == x1)
        y = (y0+y1)/2;
    else {
        assert(x1 >= x0);
        y = y0 + (y1 - y0) * (x - x0)/(x1 - x0);
    }

    return y;
}
void gridpp::util::not_implemented_error() {
    gridpp::util::error("Not implemented");
}
void gridpp::util::check_equal_size(const vec& v1, const vec& v2) {
    if(v1.size() != v2.size()) {
        gridpp::util::error("Vectors are not of the same size");
    }
}
