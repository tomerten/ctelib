#ifndef PHYSICS_H
#define PHYSICS_H
#include <map>
#include <string>
#include <vector>
namespace PHYSICS {
std::vector<double> vectorMultiply(std::vector<double> x,
                                   std::vector<double> y);

std::vector<double> vectorAdd(std::vector<double> x, std::vector<double> y);

void updateRad(std::vector<std::vector<double>> &dist,
               std::map<std::string, double> &radparam, int seed);

void updateRF(std::vector<std::vector<double>> &dist, double fmix,
              std::map<std::string, double> &tw,
              std::map<std::string, double> &longparam, std::vector<double> &hs,
              std::vector<double> &vs);

void updateBeta(std::vector<std::vector<double>> &dist,
                std::map<std::string, double> tw, double coupling, double K2L,
                double K2SL);

} // namespace PHYSICS
#endif