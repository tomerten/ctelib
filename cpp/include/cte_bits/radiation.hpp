
#ifndef RADIATION_H
#define RADIATION_H
// #include <CL/cl.h>
#include <map>
#include <string>
#include <vector>

namespace cte_radiation {

std::map<std::string, double>
radiationEquilib(std::map<std::string, double> &twiss,
                 std::map<std::string, double> &longparam);

double RadiationLossesPerTurn(std::map<std::string, double> &twiss);
void CalcRadDecayExcitation(std::map<std::string, double> &twheader,
                            std::map<std::string, double> &radparam);

struct ranFunctor {
  std::vector<double> operator()(std::vector<double> &P, std::vector<double> &D,
                                 std::vector<double> &E,
                                 std::vector<double> &R);
};
void updateRad(std::vector<std::vector<double>> &dist,
               std::map<std::string, double> &radparam, int seed);
} // namespace cte_radiation
#endif
