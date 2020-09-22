#include "gridpp.h"
#include <iostream>

using namespace gridpp;

namespace {
    vec2 get_curve(const gridpp::Grid& rgrid, const vec3& rvalues, const vec3& n_at_r, int y, int x, float outer_radius) {
        vec2 curve;
        vec distances;
        vec2 rlats = rgrid.get_lats();
        vec2 rlons = rgrid.get_lons();
        int Y = rgrid.size()[0];
        int X = rgrid.size()[1];
        int T = rvalues.size();
        ivec2 outer = rgrid.get_neighbours_with_distance(rlats[y][x], rlons[y][x], outer_radius, distances);

        // add the current gridpoint, since the neighbour search omits this point
        ivec temp(2);
        temp[0] = y;
        temp[1] = x;
        outer.push_back(temp);
        distances.push_back(0);

        // Create a curve of observed and radar values based on the neighbourhood points
        int J = outer.size();
        vec x_curve;
        vec y_curve;
        x_curve.reserve(J);
        y_curve.reserve(J);
        float xtotal = 0;
        float ytotal = 0;
        float xtotal0 = 0;
        float ytotal0 = 0;
        int num_nonzeros = 0;
        assert(outer.size() == distances.size());
        for(int t = 0; t < T; t++) {
            for(int j = 0; j < J; j++) {
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                // if(gridpp::is_valid(n_at_r[outer[j]]) && rvalues[outer[j]] > 0 && n_at_r[outer[j]] > 0) {
                // Check that we have a point with enough observations
                float xx = rvalues[t][y_outer][x_outer];
                float yy = n_at_r[t][y_outer][x_outer];
                if(gridpp::is_valid(xx) && gridpp::is_valid(yy)) {
                    x_curve.push_back(xx);
                    y_curve.push_back(yy);
                    xtotal += xx;
                    ytotal += yy;
                    if(xx > 0 && yy > 0) {
                        xtotal0 += xx;
                        ytotal0 += yy;
                        num_nonzeros++;
                    }
                }
            }
        }
        x_curve.push_back(0);
        y_curve.push_back(0);
        std::sort(x_curve.begin(), x_curve.end());
        std::sort(y_curve.begin(), y_curve.end());
        curve.push_back(x_curve);
        curve.push_back(y_curve);
        return curve;
    }
}

vec2 gridpp::correction(const Grid& rgrid,
        const vec3& rvalues,
        const Points& npoints,
        const vec2& nvalues,
        const vec2& apply_values,
        float mean_radius,
        float outer_radius,
        float inner_radius,
        int min_num,
        int max_num,
        gridpp::CorrectionType type,
        ivec2& count) {
    int Y = rgrid.size()[0];
    int X = rgrid.size()[1];
    int T = rvalues.size();

    vec2 output = gridpp::init_vec2(Y, X, 0);

    // Store the number of observations used to correct each gridpoint
    count = gridpp::init_ivec2(Y, X, -1);

    vec2 rlats = rgrid.get_lats();
    vec2 rlons = rgrid.get_lons();

    // Compute the average observed at each radar point
    vec3 n_at_r(T);
    for(int t = 0; t < T; t++) {
        n_at_r[t] = gridpp::gridding(rgrid, npoints, nvalues[t], mean_radius, 1, gridpp::Mean);
    }
    // vec radii(npoints.size(), outer_radius);
    // vec2 zeros = gridpp::init_vec2(Y, X, 0);
    // vec2 has_obs = gridpp::fill(rgrid, zeros, npoints, radii, 1, false);

    int num_corrected = 0;

    // Count the number of times (search radii) we compute calibration coefficients
    int iterations = 0;

    for(int y = 0; y < Y; y++) {
        if(y % 100 == 0)
            std::cout << y << "/" << Y << std::endl;
        for(int x = 0; x < X; x++) {
            // if(has_obs[y][x] == 0)
            //     // No point trying to create calibration, because there are no nearby observations
            //     continue;

            if(count[y][x] > -1)
                // This mean we have already processed this gridpoint (because it was included in
                // the search radius of another grid point that was processed earlier)
                continue;
            iterations++;

            // Retrive a list of nearby points and distances
            vec distances;
            ivec2 outer = rgrid.get_neighbours_with_distance(rlats[y][x], rlons[y][x], outer_radius, distances);

            // add the current gridpoint, since the neighbour search omits this point
            ivec temp(2);
            temp[0] = y;
            temp[1] = x;
            outer.push_back(temp);
            distances.push_back(0);

            // Create a curve of observed and radar values based on the neighbourhood points
            int J = outer.size();
            vec x_curve;
            vec y_curve;
            x_curve.reserve(J);
            y_curve.reserve(J);
            float xtotal = 0;
            float ytotal = 0;
            float xtotal0 = 0;
            float ytotal0 = 0;
            int num_nonzeros = 0;
            assert(outer.size() == distances.size());
            for(int t = 0; t < T; t++) {
                for(int j = 0; j < J; j++) {
                    int y_outer = outer[j][0];
                    int x_outer = outer[j][1];
                    // if(gridpp::is_valid(n_at_r[outer[j]]) && rvalues[outer[j]] > 0 && n_at_r[outer[j]] > 0) {
                    // Check that we have a point with enough observations
                    float xx = rvalues[t][y_outer][x_outer];
                    float yy = n_at_r[t][y_outer][x_outer];
                    if(gridpp::is_valid(xx) && gridpp::is_valid(yy)) {
                        x_curve.push_back(xx);
                        y_curve.push_back(yy);
                        xtotal += xx;
                        ytotal += yy;
                        if(xx > 0 && yy > 0) {
                            xtotal0 += xx;
                            ytotal0 += yy;
                            num_nonzeros++;
                        }
                    }
                }
            }
            x_curve.push_back(0);
            y_curve.push_back(0);
            std::sort(x_curve.begin(), x_curve.end());
            std::sort(y_curve.begin(), y_curve.end());

            // Find the threshold where should truncate precip to 0
            float x_threshold = 0;
            for(int q = 0; q < x_curve.size(); q++) {
                if(y_curve[q] > 0) {
                    // std::cout << q << " " << x_curve[q] << " " << y_curve[q] << std::endl;
                    // std::cout << xtotal / x_curve.size() << " " << ytotal / y_curve.size() << std::endl;
                    x_threshold = x_curve[q-1];
                    break;
                }
            }
            int S = x_curve.size();

            // TODO: Deal with case where original is 0. We shouldn't generate precip then, or?

            // Remove the extremes of the curve
            if(S > 20) {
                x_curve = vec(x_curve.begin(), x_curve.end() - S / 10);
                y_curve = vec(y_curve.begin(), y_curve.end() - S / 10);
            }
            // List of values in neighbourhood that we want to adjust
            vec curr_values(J);
            for(int j = 0; j < J; j++) {
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                curr_values[j] = apply_values[y_outer][x_outer];
            }

            // Compute adjusted values
            vec curr_adjusted_values;
            if(type == gridpp::Qq) {
                vec2 curve = gridpp::quantile_mapping_curve(y_curve, x_curve);
                curr_adjusted_values = gridpp::apply_curve(curr_values, curve, gridpp::Zero, gridpp::Zero);
            }
            else {
                if(type == gridpp::Multiplicative) {
                    curr_adjusted_values.resize(J);
                    float ratio = 1;
                    if(num_nonzeros >= min_num)
                        ratio = ytotal / xtotal;
                    for(int j = 0; j < J; j++)
                        curr_adjusted_values[j] = curr_values[j] * ratio;
                }
                else {
                    curr_adjusted_values.resize(J);
                    float ratio = 1;
                    if(num_nonzeros >= min_num)
                        ratio = ytotal0 / xtotal0;
                    for(int j = 0; j < J; j++) {
                        if(curr_values[j] < x_threshold)
                            curr_adjusted_values[j] = 0;
                        else
                            curr_adjusted_values[j] = curr_values[j] * ratio;
                    }
                }

            }

            // Apply to values inside the circle
            for(int j = 0; j < J; j++) {
                assert(j < outer.size());
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                int curr_count = count[y_outer][x_outer];
                float curr_r = apply_values[y_outer][x_outer];
                // Don't update this gridpoint if there is less info now
                if(curr_count >= S)
                    continue;
                if(distances[j] < inner_radius) {
                    count[y_outer][x_outer] = S;
                    if(S >= min_num) {
                        output[y_outer][x_outer] = curr_adjusted_values[j];
                        num_corrected++;
                    }
                    else {
                        output[y_outer][x_outer] = curr_r;
                    }
                }
            }
        }
    }
    std::cout << "Corrected: " << num_corrected << " iterations: " << iterations << std::endl;
    return output;
}
vec2 gridpp::correction(const Grid& rgrid,
        const vec2& rvalues,
        const Points& npoints,
        const vec& nvalues,
        float mean_radius,
        float outer_radius,
        float inner_radius,
        int min_num,
        int max_num,
        gridpp::CorrectionType type,
        ivec2& count) {
    vec3 rvalues3(1);
    rvalues3[0] = rvalues;
    vec2 nvalues2(1);
    nvalues2[0] = nvalues;
    vec2 results = correction(rgrid, rvalues3, npoints, nvalues2, rvalues, mean_radius, outer_radius, inner_radius, min_num, max_num, type, count);
    return results;
    /*
    int Y = rgrid.size()[0];
    int X = rgrid.size()[1];

    vec2 output = gridpp::init_vec2(Y, X, gridpp::MV);
    count = gridpp::init_vec2(Y, X, -1);
    vec2 rlats = rgrid.get_lats();
    vec2 rlons = rgrid.get_lons();

    // Compute the average observed at each radar point
    vec2 n_at_r = gridpp::gridding(rgrid, npoints, nvalues, mean_radius, 1, gridpp::Mean);

    int num_corrected = 0;
    int iterations = 0;

    for(int y = 0; y < Y; y++) {
        std::cout << y << "/" << Y << std::endl;
        for(int x = 0; x < X; x++) {
            if(x == 50) {
                std::cout << y << " " << x << " " << count[y][x] << std::endl;
            }
            if(count[y][x] > 0)
                continue;
            iterations++;

            // Retrive a list of nearby points and distances
            vec distances;
            ivec2 outer = rgrid.get_neighbours_with_distance(rlats[y][x], rlons[y][x], outer_radius, distances);

            // add the current gridpoint, since the neighbour search omits this point
            ivec temp(2);
            temp[0] = y;
            temp[1] = x;
            outer.push_back(temp);
            distances.push_back(0);

            // Create a curve of observed and radar values based on the neighbourhood points
            int J = outer.size();
            vec x_curve;
            vec y_curve;
            x_curve.reserve(J);
            y_curve.reserve(J);
            float xtotal = 0;
            float ytotal = 0;
            float xtotal0 = 0;
            float ytotal0 = 0;
            int num_nonzeros = 0;
            assert(outer.size() == distances.size());
            for(int j = 0; j < J; j++) {
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                // if(gridpp::is_valid(n_at_r[outer[j]]) && rvalues[outer[j]] > 0 && n_at_r[outer[j]] > 0) {
                // Check that we have a point with enough observations
                assert(y_outer < rvalues.size());
                assert(y_outer < n_at_r.size());
                assert(x_outer < rvalues[y_outer].size());
                assert(x_outer < n_at_r[y_outer].size());
                float xx = rvalues[y_outer][x_outer];
                float yy = n_at_r[y_outer][x_outer];
                if(gridpp::is_valid(xx) && gridpp::is_valid(yy)) {
                    x_curve.push_back(xx);
                    y_curve.push_back(yy);
                    xtotal += xx;
                    ytotal += yy;
                    if(xx > 0 && yy > 0) {
                        xtotal0 += xx;
                        ytotal0 += yy;
                        num_nonzeros++;
                    }
                }
            }
            x_curve.push_back(0);
            y_curve.push_back(0);
            std::sort(x_curve.begin(), x_curve.end());
            std::sort(y_curve.begin(), y_curve.end());

            // Find the threshold where should truncate precip to 0
            float x_threshold = 0;
            for(int q = 0; q < x_curve.size(); q++) {
                if(y_curve[q] > 0) {
                    // std::cout << q << " " << x_curve[q] << " " << y_curve[q] << std::endl;
                    // std::cout << xtotal / x_curve.size() << " " << ytotal / y_curve.size() << std::endl;
                    x_threshold = x_curve[q-1];
                    break;
                }
            }
            if(x == 50) {
                std::cout << y << " " << x << " " << x_threshold << std::endl;
            }
            int S = x_curve.size();

            // Remove the extremes of the curve
            if(S > 20) {
                x_curve = vec(x_curve.begin(), x_curve.end() - S / 10);
                y_curve = vec(y_curve.begin(), y_curve.end() - S / 10);
            }
            if(J == 0) {
                output[y][x] = rvalues[y][x];
                count[y][x] = 0;
            }
            vec applied;
            vec curr_rvalues(J);
            for(int j = 0; j < J; j++) {
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                curr_rvalues[j] = rvalues[y_outer][x_outer];
            }
            if(type == gridpp::Qq) {
                vec2 curve = gridpp::quantile_mapping_curve(y_curve, x_curve);
                applied = gridpp::apply_curve(curr_rvalues, curve, gridpp::Zero, gridpp::Zero);
            }
            else {
                if(type == gridpp::Multiplicative) {
                    applied.resize(J);
                    float ratio = 1;
                    if(num_nonzeros >= min_num)
                        ratio = ytotal / xtotal;
                    std::cout << y << " " << " " << ratio << std::endl;
                    assert(ratio >= 0);
                    assert(ratio <= 100);
                    for(int j = 0; j < J; j++)
                        applied[j] = curr_rvalues[j] * ratio;
                }
                else {
                    applied.resize(J);
                    float ratio = 1;
                    if(num_nonzeros >= min_num)
                        ratio = ytotal0 / xtotal0;
                    assert(ratio >= 0);
                    assert(ratio <= 100);
                    for(int j = 0; j < J; j++) {
                        if(curr_rvalues[j] < x_threshold)
                            applied[j] = 0;
                        else
                            applied[j] = curr_rvalues[j] * ratio;
                    }
                }

            }

            for(int j = 0; j < J; j++) {
                assert(j < outer.size());
                int y_outer = outer[j][0];
                int x_outer = outer[j][1];
                float curr_count = count[y_outer][x_outer];
                float curr_r = rvalues[y_outer][x_outer];
                if(curr_count >= num_nonzeros)
                    continue;
                if(distances[j] < inner_radius) {
                    curr_count = num_nonzeros;
                    count[y_outer][x_outer] = curr_count;
                    // curr_count = ytotal / x_curve.size() - xtotal / x_curve.size();
                    if(num_nonzeros >= min_num) {
                        if(1 || curr_r > 0.1) {
                            output[y_outer][x_outer] = applied[j];
                            num_corrected++;
                        }
                        else {
                            output[y_outer][x_outer] = curr_r;
                        }
                    }
                    else {
                        output[y_outer][x_outer] = curr_r;
                        count[y_outer][x_outer] = 1;
                    }
                }
            }
            // std::cout << "Size: " << x_curve.size() << " " << xtotal / x_curve.size() << " " << ytotal / y.size() << std::endl;
        }
    }
    std::cout << "Corrected: " << num_corrected << " iterations: " << iterations << std::endl;
    return output;
    */
}
