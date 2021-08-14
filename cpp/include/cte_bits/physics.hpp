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

void updateIBS(std::vector<std::vector<double>> &dist,
               std::map<std::string, double> &lparam,
               std::map<std::string, double> &bparam,
               std::map<std::string, double> &twheader, std::vector<double> &h,
               std::vector<double> &ibsCoeffLast,
               std::vector<double> &ibsHistCoeffLast);

} // namespace PHYSICS
#endif