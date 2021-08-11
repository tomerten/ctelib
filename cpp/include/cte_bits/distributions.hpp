#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
// #include <CL/cl.h>
#include <map>
#include <string>
#include <vector>

namespace cte_distributions {
std::vector<double> BiGaussian4D(const double betax, const double ex,
                                 const double betay, const double ey, int seed);

std::vector<double>
BiGaussian6DLongMatched(double betax, double ex, double betay, double ey,
                        std::vector<double> &h, std::vector<double> &v,
                        std::map<std::string, double> &twheader,
                        std::map<std::string, double> &longparam, int seed);

std::vector<std::vector<double>>
GenerateDistributionMatched(int nMacro, double betax, double ex, double betay,
                            double ey, std::vector<double> &h,
                            std::vector<double> &v,
                            std::map<std::string, double> &twheader,
                            std::map<std::string, double> &longparam, int seed);

} // namespace cte_distributions

#endif