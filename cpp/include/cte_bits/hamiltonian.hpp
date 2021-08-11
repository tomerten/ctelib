#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
// #include <CL/cl.h>
#include <map>
#include <string>
#include <vector>

namespace cte_hamiltonian {
double tcoeff(std::map<std::string, double> &twissheaderL,
              double baseHarmonicNumber);

double pcoeff(std::map<std::string, double> &twissheaderL, double voltage);

double Hamiltonian(std::map<std::string, double> &twheader,
                   std::map<std::string, double> &longparam,
                   std::vector<double> &h, std::vector<double> &v,
                   double tcoeff, double t, double delta);
} // namespace cte_hamiltonian

#endif