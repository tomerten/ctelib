#ifndef CTESYNCH_H
#define CTESYNCH_H
#include <algorithm>
#include <map>
#include <math.h>
#include <string>
#include <vector>

namespace CTESYNCH {
template <class T> struct synchronousPhaseFunctor {
  synchronousPhaseFunctor(T const &target, std::vector<double> &voltages,
                          std::vector<double> &harmonicNumbers, double charge)
      : U0(target), volts(voltages), hs(harmonicNumbers), ch(charge) {}
  std::tuple<double, double> operator()(T const &phi);

private:
  T U0;
  std::vector<double> volts;
  std::vector<double> hs;
  double ch;
};

template <class T>
T synchronousPhaseFunctorRoot(T x, std::vector<double> &voltages,
                              std::vector<double> &harmnumbers, double charge,
                              T guess, T min, T max);

double VoltageRf(double phi, std::vector<double> &volts,
                 std::vector<double> &hs);

double VoltageRfPrime(double phi, double charge, std::vector<double> &volts,
                      std::vector<double> &hs);

double SynchrotronTune(std::map<std::string, double> &twheader,
                       std::map<std::string, double> &longparams,
                       std::vector<double> &volts, std::vector<double> &hs);
} // namespace CTESYNCH

#endif