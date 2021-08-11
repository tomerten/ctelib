#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <cstdlib>
#include <ibs>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

namespace CTESYNCH {
/*
 *****************************************************************************
 *****************************************************************************
 * FUNCTOR TO RETURN VOLTAGE * CHARGE - U0 AND DERIVATIVE AS TUPLE
 * THIS IS USED AS INPUT FOR BOOST NEWTON RAPHSON ROOT SEARCH.
 *****************************************************************************
 * Authors:
 *  - Tom Mertens
 *
 * History:
 *  - 10/08/2021 : updated version from original ste code
 *****************************************************************************
 * Arguments:
 * ----------
 *  - T const &target
 *      Target value for the voltage (RF compensate for U0 so U0 in eV)
 *  - std::vector<double> &voltages
 *      Voltages of RF systems IMPORTANT: first element is for main RF
 *  - std::vector<double> &harmonicNumbers
 *      Harmonic numbers  IMPORTANT: first element is for main RF
 *  - double charge
 *      Particle charge
 *  - T const &phi
 *      phase in rad
 ******************************************************************************
 ******************************************************************************
 */
template <class T> struct synchronousPhaseFunctor {
  synchronousPhaseFunctor(T const &target, std::vector<double> &voltages,
                          std::vector<double> &harmonicNumbers, double charge)
      : U0(target), volts(voltages), hs(harmonicNumbers), ch(charge) {}
  std::tuple<double, double> operator()(T const &phi) {

    // init
    T vrf = ch * volts[0] * sin(phi);
    T dvrf = ch * volts[0] * cos(phi);

    // add the rest taking harmonic numbers into account
    for (int i = 1; i < hs.size(); i++) {
      vrf += ch * volts[i] * sin((hs[i] / hs[0]) * phi);
      dvrf += ch * volts[i] * (hs[i] / hs[0]) * cos((hs[i] / hs[0]) * phi);
    }

    std::tuple<double, double> out = {vrf - U0, dvrf};
    return out;
  }

private:
  T U0;
  std::vector<double> volts;
  std::vector<double> hs;
  double ch;
};

/*
================================================================================
================================================================================
BOOST NEWTON RAPHSON ROOT SEARCH FOR SYNCHRONOUS PHASE.
================================================================================

================================================================================
Arguments:

//phi is in rad
================================================================================
================================================================================
*/
template <class T>
T synchronousPhaseFunctorRoot(T x, std::vector<double> &voltages,
                              std::vector<double> &harmnumbers, double charge,
                              T guess, T min, T max) {
  // return cube root of x using 1st derivative and Newton_Raphson.
  using namespace boost::math::tools;

  const int digits =
      std::numeric_limits<T>::digits; // Maximum possible binary digits accuracy
                                      // for type T.
  int get_digits = static_cast<int>(
      digits * 0.6); // Accuracy doubles with each step, so stop when we have
                     // just over half the digits correct.
  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  T result = newton_raphson_iterate(
      synchronousPhaseFunctor<T>(x, voltages, harmnumbers, charge), guess, min,
      max, get_digits, it);
  return result;
};

/********************************************************************************
 ********************************************************************************
 * CALCULATE TOTAL RF VOLTAGE FOR GIVEN PHASE (IN RAD)
 *
 ********************************************************************************
 */
double VoltageRf(double phi, std::vector<double> &volts,
                 std::vector<double> &hs) {
  double vrf = volts[0] * sin(phi);

  for (int i = 1; i < hs.size(); i++) {
    vrf += volts[i] * sin((hs[i] / hs[0]) * phi);
  }

  return vrf;
};

double VoltageRfPrime(double phi, double charge, std::vector<double> &volts,
                      std::vector<double> &hs) {
  // init
  double vrf = volts[0] * cos(phi);

  // add other rfs
  for (int i = 1; i < volts.size(); i++) {
    vrf += volts[i] * (hs[i] / hs[0]) * cos((hs[i] / hs[0]) * phi);
  }

  // V -> eV
  vrf *= charge;
  return vrf;
}

double SynchrotronTune(std::map<std::string, double> &twheader,
                       std::map<std::string, double> &longparams,
                       std::vector<double> &volts, std::vector<double> &hs) {

  double charge = twheader["CHARGE"];
  double phis = longparams["phis"];
  double pc = twheader["PC"];

  return sqrt(hs[0] * twheader["eta"] *
              fabs(charge * VoltageRfPrime(phis, charge, volts, hs)) /
              (2 * pi * pc * 1e9));
}
} // namespace CTESYNCH