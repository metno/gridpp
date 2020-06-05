#include "gridpp.h"
#include <iostream>

vec gridpp::correction(const Points& rpoints, const vec& rvalues, const Points& npoints, const vec& nvalues, float mean_radius, float outer_radius, float inner_radius, int min_num, int max_num, std::string type, vec& count) {
    int s = rpoints.size();
    int itype = 0;
    if(type == "qq")
        itype = 0;
    else if(type == "mult")
        itype = 1;
    else
        itype = 2;
    vec output(s, 0);
    count.clear();
    count.resize(s, -1);
    vec n_at_r(s, gridpp::MV);
    double s_time = gridpp::util::clock();
    int num_points = 0;
    vec rlats = rpoints.get_lats();
    vec rlons = rpoints.get_lons();
    for(int i = 0; i < s; i++) {
        ivec I = npoints.get_neighbours(rlats[i], rlons[i], mean_radius);
        float total = 0;
        int c = 0;
        if(I.size() > 0) {
            for(int j = 0; j < I.size(); j++) {
                assert(nvalues.size() > I[j]);
                total += nvalues[I[j]];
                c++;
            }
            if(c > 0) {
                n_at_r[i] = total / c;
                num_points++;
            }
        }
    }
    double e_time = gridpp::util::clock();
    std::cout << e_time - s_time << std::endl;
    std::cout << num_points << std::endl;
    int num_corrected = 0;
    int iterations = 0;

    for(int i = 0; i < s; i++) {
        if(i % int(s / 20) == 0)
            std::cout << i << "/" << s << std::endl;
        if(count[i] > 0)
            continue;
        iterations++;
        vec distances;
        ivec outer = rpoints.get_neighbours_with_distance(rlats[i], rlons[i], outer_radius, distances);
        outer.push_back(i);
        distances.push_back(0);

        int J = outer.size();
        vec x;
        vec y;
        x.reserve(J);
        y.reserve(J);
        float xtotal = 0;
        float ytotal = 0;
        float xtotal0 = 0;
        float ytotal0 = 0;
        int num_nonzeros = 0;
        assert(outer.size() == distances.size());
        for(int j = 0; j < J; j++) {
            // if(gridpp::util::is_valid(n_at_r[outer[j]]) && rvalues[outer[j]] > 0 && n_at_r[outer[j]] > 0) {
            if(gridpp::util::is_valid(n_at_r[outer[j]])) {
                x.push_back(rvalues[outer[j]]);
                y.push_back(n_at_r[outer[j]]);
                xtotal += rvalues[outer[j]];
                ytotal += n_at_r[outer[j]];
                if(rvalues[outer[j]] > 0 && n_at_r[outer[j]] > 0) {
                    xtotal0 += rvalues[outer[j]];
                    ytotal0 += n_at_r[outer[j]];
                    num_nonzeros++;
                }
            }
        }
        x.push_back(0);
        y.push_back(0);
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        float x_threshold = 0;
        for(int q = 0; q < x.size(); q++) {
            if(y[q] > 0) {
                // std::cout << q << " " << x[q] << " " << y[q] << std::endl;
                // std::cout << xtotal / x.size() << " " << ytotal / y.size() << std::endl;
                x_threshold = x[q-1];
                break;
            }
        }
        int X = x.size();
        if(X > 20) {
            x = vec(x.begin(), x.end() - X / 10);
            y = vec(y.begin(), y.end() - X / 10);
        }
        if(J == 0) {
            output[i] = rvalues[i];
            count[i] = 0;
        }
        for(int j = 0; j < J; j++) {
            int index = outer[j];
            if(count[index] >= num_nonzeros)
                continue;
            if(distances[j] < inner_radius) {
                count[index] = num_nonzeros;
                // count[index] = ytotal / x.size() - xtotal / x.size();
                if(num_nonzeros >= min_num) {
                    if(1 || rvalues[index] > 0.1) {
                        if(itype == 0)
                            output[index] = gridpp::quantile_mapping(rvalues[index], x, y, gridpp::Zero);
                        else if (itype == 1)
                            output[index] = rvalues[index] * ytotal / xtotal;
                        else {
                            if(rvalues[index] < x_threshold)
                                output[index] = 0;
                            else if(xtotal0 > 0){
                                output[index] = rvalues[index] * ytotal0 / xtotal0;
                            }
                            else {
                                output[index] = rvalues[index];
                            }
                        }
                        num_corrected++;
                    }
                    else {
                        output[index] = rvalues[index];
                    }
                }
                else {
                    output[index] = rvalues[index];
                    count[index] = 1;
                }
            }
        }
        // std::cout << "Size: " << x.size() << " " << xtotal / x.size() << " " << ytotal / y.size() << std::endl;
    }
    std::cout << "Corrected: " << num_corrected << " iterations: " << iterations << std::endl;
    return output;
}
