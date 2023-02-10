#include "gridpp.h"
#include <math.h>
#include <iostream>

using namespace gridpp;

// Notes:
// Bilinear should have roughly the same performance as nearest. THe exception is for 3D arrays
// where the number of timesteps is large. In this case, nearest processes time in the outer loop to
// prevent inefficient memory access. This is too complicated to implement here at the moment.

namespace {
    // Bilinear interpolation for a given point
    float calc(const Grid& grid, const vec2& iInputLats, const vec2& iInputLons, const vec2& ivalues, float lat, float lon);
    vec calc(const Grid& grid, const vec2& iInputLats, const vec2& iInputLons, const vec3& ivalues, float lat, float lon);
    // Bilinear interpolation based on 4 surrounding points with coordinates (x0,y0), (x1,y1), etc
    // and values v0, v1, etc
    float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3);
    // Compute s,t when the points form a parallelogram
    bool calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s);

    // Check that value is within -tol and 1+tol
    bool is_within_range(float value);
    // Compute s,t when the points do not form a parallelogram
    bool calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s);
}
vec2 gridpp::bilinear(const Grid& igrid, const Grid& ogrid, const vec2& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output = gridpp::init_vec2(nLat, nLon, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    // Algorithm from here:
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            // Use the four points surrounding the lookup point. Use the nearest neighbour
            // and then figure out what side of this point the other there points are.
            float lat = iOutputLats[i][j];
            float lon = iOutputLons[i][j];
            output[i][j] = ::calc(igrid, iInputLats, iInputLons, ivalues, lat, lon);
        }
    }
    return output;
}
vec3 gridpp::bilinear(const Grid& igrid, const Grid& ogrid, const vec3& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
       throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nTime = ivalues.size();
    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec3 output = gridpp::init_vec3(nTime, nLat, nLon, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    // To reuse nearest neighbour information across time, we can't call calc on each timestep
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            float lat = iOutputLats[i][j];
            float lon = iOutputLons[i][j];
            vec temp = ::calc(igrid, iInputLats, iInputLons, ivalues, lat, lon);
            for(int t = 0; t < nTime; t++)
                output[t][i][j] = temp[t];
        }
    }
    return output;
}

vec gridpp::bilinear(const Grid& igrid, const Points& opoints, const vec2& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    // Algorithm from here:
    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        // Use the four points surrounding the lookup point. Use the nearest neighbour
        // and then figure out what side of this point the other there points are.
        float lat = iOutputLats[i];
        float lon = iOutputLons[i];
        output[i] = ::calc(igrid, iInputLats, iInputLons, ivalues, lat, lon);
    }
    return output;
}
vec2 gridpp::bilinear(const Grid& igrid, const Points& opoints, const vec3& ivalues) {
    if(!gridpp::compatible_size(igrid, ivalues))
       throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nTime = ivalues.size();
    int nPoints = iOutputLats.size();

    vec2 output = gridpp::init_vec2(nTime, nPoints, gridpp::MV);
    if(igrid.size()[0] == 0 || igrid.size()[1] == 0)
        return output;

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        float lat = iOutputLats[i];
        float lon = iOutputLons[i];
        vec temp = ::calc(igrid, iInputLats, iInputLons, ivalues, lat, lon);
        for(int t = 0; t < nTime; t++) {
            output[t][i] = temp[t];
        }
    }
    return output;
}

namespace {
    bool calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s) {
        // std::cout << "Method 3: Parallelogram" << std::endl;
        float X31 = X3 - X1;
        float X21 = X2 - X1;
        float Y21 = Y2 - Y1;
        float Y31 = Y3 - Y1;

        float A = X21;
        float B = X31;
        float C = Y21;
        float D = Y31;
        assert(A * D != B * C);
        float det = 1 / (A * D - B * C);
        s = det * ((x - X1) * (D) + (y - Y1) * (-B));
        t = det * ((x - X1) * (-C) + (y - Y1) * (A));
        return true;
    }
    bool is_within_range(float value) {
        float tol = 0.01;
        return value >= -tol && value < 1 + tol;
    }

    bool calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s) {
        // std::cout << "Method general" << std::endl;
        double a = -x0 + x2;
        double b = -x0 + x1;
        double c = x0 - x1 - x2 + x3;
        double d = x - x0;
        double e = -y0 + y2;
        double f = -y0 + y1;
        double g = y0 - y1 - y2 + y3;
        double h = y - y0;
        double alpha=gridpp::MV, beta=gridpp::MV;
        double Y1 = y1;
        double Y2 = y3;
        double Y3 = y0;
        double Y4 = y2;
        double X1 = x1;
        double X2 = x3;
        double X3 = x0;
        double X4 = x2;
        double X31 = X3 - X1;
        double X21 = X2 - X1;
        double Y42 = Y4 - Y2;
        double Y21 = Y2 - Y1;
        double Y31 = Y3 - Y1;
        double Y43 = Y4 - Y3;
        double X42 = X4 - X2;
        double X43 = X4 - X3;
        double tol = 0.001;
        /*
        std::cout << "a " << a << std::endl;
        std::cout << "b " << b << std::endl;
        std::cout << "c " << c << std::endl;
        std::cout << "d " << d << std::endl;
        std::cout << "e " << e << std::endl;
        std::cout << "f " << f << std::endl;
        std::cout << "g " << g << std::endl;
        std::cout << "h " << h << std::endl;
  */
        if(2 * c * e - 2 * a * g != 0 && 2 * c * f - 2 * b * g != 0) {
            // Quadratic equation has two solutions
            alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
            beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
            /*
            std::cout << "alpha " << alpha << " beta " << beta << std::endl;
            std::cout << "sqrt " << sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)) << std::endl;
            std::cout << "first part " << b * e - a * f + d * g - c * h << std::endl;
            std::cout << "second part " << b * e - a * f - d * g + c * h << std::endl;
            */
            if(!is_within_range(alpha))
                alpha = -(b * e - a * f + d * g - c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
            if(!is_within_range(beta))
                beta  = (b * e - a * f - d * g + c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
            if(!is_within_range(alpha) || !is_within_range(beta)) {
            }
            /*
            std::cout << std::setprecision(10) << "alpha " << alpha << " beta " << beta << std::endl;
            std::cout << std::setprecision(10) << x << std::endl;
            std::cout << std::setprecision(10) << y << std::endl;
            std::cout << std::setprecision(10) << x0 << std::endl;
            std::cout << std::setprecision(10) << x1 << std::endl;
            std::cout << std::setprecision(10) << x2 << std::endl;
            std::cout << std::setprecision(10) << x3 << std::endl;
            std::cout << std::setprecision(10) << y0 << std::endl;
            std::cout << std::setprecision(10) << y1 << std::endl;
            std::cout << std::setprecision(10) << y2 << std::endl;
            std::cout << std::setprecision(10) << y3 << std::endl;
  */
            assert(is_within_range(alpha));
            assert(is_within_range(beta));
        }
        else if(2 * c * f - 2 * b * g == 0) {
            // std::cout << "Need to diagnose t" << std::endl;
            alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
            if(!is_within_range(alpha))
                alpha = -(b * e - a * f + d * g - c * h - sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
            float s = alpha;
            float t;
            if(Y3 + Y43 * s - Y1 - Y21 * s == 0) {
                assert(X3 + X43 * s - X1 - X21 * s != 0);
                t = (x - X1 - X21 * s) / (X3 + X43 * s - X1 - X21 * s);
            }
            else {
                t = (y - Y1 - Y21 * s) / (Y3 + Y43 * s - Y1 - Y21 * s);
            }
            beta = 1 - t;
            assert(is_within_range(alpha));
            assert(is_within_range(beta));
        }
        else if(2 * c * e - 2 * a * g == 0) {
            // std::cout << "Need to diagnose s" << std::endl;
            beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
            if(!is_within_range(beta))
                beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
            float t =  1 - beta;
            float s;
            if(Y2 + Y42 * t - Y1 - Y31 * t == 0) {
                assert(X2 + X42 * t - X1 - X31 * t != 0);
                s = (x - X1 - X31 * t) / (X2 + X42 * t - X1 - X31 * t);
            }
            else {
                s = (y - Y1 - Y31 * t) / (Y2 + Y42 * t - Y1 - Y31 * t);
            }
            alpha = s;
            assert(is_within_range(alpha));
            assert(is_within_range(beta));
        }
        s = alpha;
        t = 1 - beta;
        return true;
    }
    float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3) {
       float Y1 = y1;
       float Y2 = y3;
       float Y3 = y0;
       float Y4 = y2;
       float X1 = x1;
       float X2 = x3;
       float X3 = x0;
       float X4 = x2;
       float P1 = v1;
       float P2 = v3;
       float P3 = v0;
       float P4 = v2;

       // General method based on: https://stackoverflow.com/questions/23920976/bilinear-interpolation-with-non-aligned-input-points
       // Parallelogram method based on: http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/InterpolationIrregularBilinear.pdf

       float s = gridpp::MV, t = gridpp::MV;
       bool rectangularGrid = (X1 == X3 && X2 == X4 && Y1 == Y2 && Y3 == Y4);
       // TODO: Why are the tolerances so high?
       bool verticalParallel = fabs((X3 - X1)*(Y4 - Y2) - (X4 - X2)*(Y3 - Y1)) <= 1e-4;
       bool horizontalParallel = fabs((X2 - X1)*(Y4 - Y3) - (X4 - X3)*(Y2 - Y1)) <= 1e-4;

       // std::cout << rectangularGrid << " " << fabs((X3 - X1)*(Y4 - Y2) - (X4 - X2)*(Y3 - Y1)) << " " << fabs((X2 - X1)*(Y4 - Y3) - (X4 - X3)*(Y2 - Y1)) << std::endl;

       // Compute s and t
       if(verticalParallel && horizontalParallel)
          calcParallelogram(x, y, X1, X2, X3, X4, Y1, Y2, Y3, Y4, t, s);
       else
          calcGeneral(x, y, x0, x1, x2, x3, y0, y1, y2, y3, t, s);

       if(t >= 1 && t <= 1.15)
          t = 1;
       if(t <= 0 && t >= -0.15)
          t = 0;
       if(s >= 1 && s <= 1.15)
          s = 1;
       if(s <= 0 && s >= -0.15)
          s = 0;
       if(!(s >= 0 && s <= 1 && t >= 0 && t <= 1)) {
          std::stringstream ss;
          ss << "Problem with bilinear interpolation. Grid is rotated/distorted in a way that is not supported. s=" << s << " and t=" << t << " are outside [-0.05,1.05].";
          throw std::runtime_error(ss.str());
       }
       assert(s >= 0 && s <= 1 && t >= 0 && t <= 1);
       float value = P1 * (1 - s) * ( 1 - t) + P2 * s * (1 - t) + P3 * (1 - s) * t + P4 * s * t;

       return value;
    }

    // Pass in input lats and lons, because otherwise we lose a lot of speed
    float calc(const Grid& grid, const vec2& iInputLats, const vec2& iInputLons, const vec2& ivalues, float lat, float lon) {
        int I1 = gridpp::MV;
        int I2 = gridpp::MV;
        int J1 = gridpp::MV;
        int J2 = gridpp::MV;
        float output = gridpp::MV;
        bool inside = grid.get_box(lat, lon, I1, J1, I2, J2);
        // std::cout << "Coords: " << I1 << " " << J1 << " " << I2 << " " << J2 << " " << inside << std::endl;
        if(inside) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            float v0 = ivalues[I1][J1];
            float v1 = ivalues[I2][J1];
            float v2 = ivalues[I1][J2];
            float v3 = ivalues[I2][J2];
            if(gridpp::is_valid(v0) && gridpp::is_valid(v1) && gridpp::is_valid(v2) && gridpp::is_valid(v3)) {
                float value = ::bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                output = value;
            }
            else {
                ivec nn = grid.get_nearest_neighbour(lat, lon);
                output = ivalues[nn[0]][nn[1]];
            }
        }
        else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            ivec nn = grid.get_nearest_neighbour(lat, lon);
            output = ivalues[nn[0]][nn[1]];
        }
        return output;
    }
    vec calc(const Grid& grid, const vec2& iInputLats, const vec2& iInputLons, const vec3& ivalues, float lat, float lon) {
        int I1 = gridpp::MV;
        int I2 = gridpp::MV;
        int J1 = gridpp::MV;
        int J2 = gridpp::MV;
        int T = ivalues.size();
        vec output(T, gridpp::MV);
        bool inside = grid.get_box(lat, lon, I1, J1, I2, J2);
        // std::cout << "Coords: " << I1 << " " << J1 << " " << I2 << " " << J2 << " " << inside << std::endl;
        if(inside) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            for(int t = 0; t < T; t++) {
                float v0 = ivalues[t][I1][J1];
                float v1 = ivalues[t][I2][J1];
                float v2 = ivalues[t][I1][J2];
                float v3 = ivalues[t][I2][J2];
                if(gridpp::is_valid(v0) && gridpp::is_valid(v1) && gridpp::is_valid(v2) && gridpp::is_valid(v3)) {
                    float value = ::bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                    output[t] = value;
                }
                else {
                    ivec nn = grid.get_nearest_neighbour(lat, lon);
                    output[t] = ivalues[t][nn[0]][nn[1]];
                }
            }
        }
        else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            ivec nn = grid.get_nearest_neighbour(lat, lon);
            for(int t = 0; t < T; t++) {
                output[t] = ivalues[t][nn[0]][nn[1]];
            }
        }
        return output;
    }
}
