#include "../include/cte_bits/bunch.hpp"
#include "../include/cte_bits/distributions.hpp"
#include "../include/cte_bits/physics.hpp"
#include "../include/cte_bits/radiation.hpp"
#include "../include/cte_bits/random.hpp"
#include "../include/cte_bits/synch.hpp"
#include "../include/cte_bits/utils.hpp"
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <fstream>
#include <ibs>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace cte_long {
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

double SynchrotronTune(double omega0, double U0, double charge,
                       std::vector<double> &volts, std::vector<double> &hs,
                       double phis, double eta, double pc) {
  return sqrt(hs[0] * eta *
              fabs(charge * VoltageRfPrime(phis, charge, volts, hs)) /
              (2 * pi * pc * 1e9));
}

} // namespace cte_long

Bunch::Bunch(std::map<std::string, double> &twissheader,
             std::map<std::string, std::vector<double>> &twiss,
             std::map<std::string, double> &bparam, std::vector<double> &h,
             std::vector<double> &v) {
  twheader = twissheader;
  tw = twiss;
  bunchParam = bparam;

  setBasic();
  setLongitudinalParameters(h, v);
  setRadiationParameters();
  setDistribution(h, v);
};

void Bunch::setBasic() {
  // basic
  double gamma = twheader["GAMMA"];
  double gammatr = twheader["GAMMATR"];
  double p0 = twheader["PC"];
  double len = twheader["LENGTH"];
  double mass = twheader["MASS"];
  double charge = twheader["CHARGE"];
  // set Atomic Number in twiss header
  twheader["aatom"] = bunchParam["atomNumber"];
  twheader["sigs"] = bunchParam["sigs"];
  twheader["betar"] = BetaRelativisticFromGamma(gamma);
  twheader["trev"] = len / (twheader["betar"] * clight);
  twheader["frev"] = 1.0 / twheader["trev"];
  twheader["omega"] = 2.0 * pi * twheader["frev"];
  twheader["eta"] = eta(gamma, gammatr);
  twheader["timeratio"] = bunchParam["timeRatio"];
}

void Bunch::setLongitudinalParameters(std::vector<double> &h,
                                      std::vector<double> &v) {
  double gamma = twheader["GAMMA"];
  double gammatr = twheader["GAMMATR"];
  double angularf = twheader["trev"] * h[0] * twheader["omega"];

  // search synch phase parameters
  double search1 = angularf / (8.0 * *std::max_element(h.begin(), h.end()));
  double search2 = search1 + angularf / *std::min_element(h.begin(), h.end());
  double searchWidth = angularf / (2.0 * *std::max_element(h.begin(), h.end()));

  // bucket number determines phis
  search1 += bunchParam["bucket"] * 2.0 * pi;
  search2 += bunchParam["bucket"] * 2.0 * pi;

  // energy loss per turn
  double U0 = cte_radiation::RadiationLossesPerTurn(twheader);

  // synchronuous phases
  double phis = cte_long::synchronousPhaseFunctorRoot(
      U0, v, h, bunchParam["charge"], search1, search1 - searchWidth,
      search1 + searchWidth);
  double phis1 = cte_long::synchronousPhaseFunctorRoot(
      U0, v, h, bunchParam["charge"], search2, search2 - searchWidth,
      search2 + searchWidth);

  double phisNext = cte_long::synchronousPhaseFunctorRoot(
      U0, v, h, bunchParam["charge"], search1 + pi,
      search1 + pi - searchWidth / 2.0, search1 + pi + searchWidth / 2.0);
  double phis1Next = cte_long::synchronousPhaseFunctorRoot(
      U0, v, h, bunchParam["charge"], search2 + pi,
      search2 + pi - searchWidth / 2.0, search2 + pi + searchWidth / 2.0);

  // save
  longitudinalParameters["sigs"] = bunchParam["sigs"];
  longitudinalParameters["phis"] = phis;
  longitudinalParameters["phis1"] = phis1;
  longitudinalParameters["phisNext"] = phisNext;
  longitudinalParameters["phisNext1"] = phis1Next;

  longitudinalParameters["qs"] =
      CTESYNCH::SynchrotronTune(twheader, longitudinalParameters, v, h);
  longitudinalParameters["tauhat"] =
      abs((phisNext - phis) / (h[0] * twheader["omega"]));

  // add to twheader for easy passing around of the values (pre-existing code)
  // twheader["phis"] = phis;
  // twheader["phisNext"] = phisNext;
  // twheader["qs"] = longitudinalParameters["qs"];
  twheader["tauhat"] = longitudinalParameters["tauhat"];
  twheader["betxavg"] = twheader["LENGTH"] / (twheader["Q1"] * 2.0 * pi);
  twheader["betyavg"] = twheader["LENGTH"] / (twheader["Q2"] * 2.0 * pi);
  longitudinalParameters["sige"] =
      sigefromsigs(twheader["omega"], bunchParam["sigs"],
                   longitudinalParameters["qs"], gamma, gammatr);
  longitudinalParameters["delta"] =
      dee_to_dpp(twheader["sige"], twheader["betar"]);
}

void Bunch::setRadiationParameters() {
  // add energy loss per turn
  double U0 = cte_radiation::RadiationLossesPerTurn(twheader);
  // add the equilibrium values and radiation damping times
  radiationParameters =
      cte_radiation::radiationEquilib(twheader, longitudinalParameters);
  radiationParameters["U0"] = U0;
  cte_radiation::CalcRadDecayExcitation(twheader, radiationParameters);
};

void Bunch::setDistribution(std::vector<double> &h, std::vector<double> &v) {
  distribution = cte_distributions::GenerateDistributionMatched(
      bunchParam["nMacro"], twheader["betxavg"], bunchParam["ex"],
      twheader["betyavg"], bunchParam["ey"], h, v, twheader,
      longitudinalParameters, bunchParam["seed"]);
}

void Bunch::getEmittance() {
  std::vector<double> out, avg;
  std::vector<double> res;
  for (int i = 0; i < distribution[0].size(); i++) {
    out.push_back(0.0);
    avg.push_back(0.0);
  }

  // calculate averages
  std::for_each(distribution.begin(), distribution.end(),
                [&](std::vector<double> v) {
                  for (int i = 0; i < distribution[0].size(); i++) {
                    avg[i] += v[i] / distribution.size();
                  }
                });

  // subtract avg, square and sum
  std::for_each(
      distribution.begin(), distribution.end(), [&](std::vector<double> v) {
        for (int i = 0; i < distribution[0].size(); i++) {
          out[i] += ((v[i] - avg[i]) * (v[i] - avg[i])) / distribution.size();
        }
      });

  out[0] /= twheader["betxavg"];
  out[2] /= twheader["betyavg"];
  out[4] = sqrt(out[4]) * clight;
  out[5] = sqrt(out[5]);

  res.push_back(out[0]);
  res.push_back(out[2]);
  res.push_back(out[4]);
  res.push_back(out[5]);

  emittances.push_back(res);
}

void Bunch::getIBSGrowthRates(int model) {

  double pnumber = bunchParam["nReal"];
  std::vector<double> emit = emittances.back();
  double ex = emit[0];
  double ey = emit[1];
  double sigs = emit[2];
  double sige = emit[3];
  double aatom = twheader["aatom"];
  double r0 = ParticleRadius(twheader["CHARGE"], aatom);
  double *ibs;
  /*
    std::printf("%-30s %12.8e\n", "aatom", aatom);
    std::printf("%-30s %12.8e\n", "ex", ex);
    std::printf("%-30s %12.8e\n", "ey", ey);
    std::printf("%-30s %12.8e\n", "ro", r0);
    std::printf("%-30s %12.8e\n", "pnumber", pnumber);
    */
  // ibs growth rates update
  switch (model) {
  case 1:
    ibs = PiwinskiSmooth(pnumber, ex, ey, sigs, sige, twheader, r0);
    break;
  case 2:
    ibs = PiwinskiLattice(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 3:
    ibs =
        PiwinskiLatticeModified(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 4:
    ibs = Nagaitsev(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 5:
    ibs =
        Nagaitsevtailcut(pnumber, ex, ey, sigs, sige, twheader, tw, r0, aatom);
    break;
  case 6:
    ibs = ibsmadx(pnumber, ex, ey, sigs, sige, twheader, tw, r0, false);
    break;
  case 7:
    ibs = ibsmadxtailcut(pnumber, ex, ey, sigs, sige, twheader, tw, r0, aatom);
    break;
  case 8:
    ibs = BjorkenMtingwa2(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 9:
    ibs = BjorkenMtingwa(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 10:
    ibs = BjorkenMtingwatailcut(pnumber, ex, ey, sigs, sige, twheader, tw, r0,
                                aatom);
    break;
  case 11:
    ibs = ConteMartini(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  case 12:
    ibs = ConteMartinitailcut(pnumber, ex, ey, sigs, sige, twheader, tw, r0,
                              aatom);
    break;
  case 13:
    ibs = MadxIBS(pnumber, ex, ey, sigs, sige, twheader, tw, r0);
    break;
  }

  std::vector<double> ibsgr;
  for (int i; i < 3; i++) {
    ibsgr.push_back(ibs[i] * 2.0 * bunchParam["fracibstot"] *
                    bunchParam["nMacro"] * bunchParam["timeRatio"] /
                    bunchParam["nReal"]);
    std::printf("%12.8e\n", ibs[i]);
  }
  ibsGrowthRates.push_back(ibsgr);
}

void Bunch::getIBSCoefficients() {
  double coeffs, coeffx, coeffy;
  double alphaAverage;
  double coeffMulT;

  std::vector<double> ibsgr = ibsGrowthRates.back();
  double alfax = ibsgr[1];
  double alfay = ibsgr[2];
  double alfap = ibsgr[0];

  double dtsamp2 = 2 * longitudinalParameters["tauhat"] / bunchParam["nbins"];
  double rmsdelta = CalcRMS(distribution)[5];
  std::vector<double> emit = emittances.back();
  double sigs = emit[3];

  double rmsx = sqrt(emit[0] * twheader["betxavg"]);
  double rmsy = sqrt(emit[1] * twheader["betyavg"]);

  // debugging
  // cout << "alfap " << alfap << endl << endl;
  if (alfap > 0.0)
    coeffs = sqrt(6.0 * alfap * twheader["trev"]) * rmsdelta;
  else
    coeffs = 0.00;

  // coupling
  if (bunchParam["ibsCoupling"] == 0.0) {
    if (alfay > 0.0)
      coeffy = sqrt(6.0 * alfay * twheader["trev"]) * rmsy;
    else
      coeffy = 0.0;

    if (alfax > 0.0)
      coeffx = sqrt(6.0 * alfax * twheader["trev"]) * rmsx;
    else
      coeffx = 0.0;
  } else {
    // alphaAverage
    alphaAverage = 0.5 * (alfax + alfay);
    if (alphaAverage > 0.0) {
      coeffx = sqrt(6.0 * alphaAverage * twheader["trev"]) * rmsx;
      coeffy = sqrt(6.0 * alphaAverage * twheader["trev"]) * rmsy;
    } else {
      coeffx = 0.0;
      coeffy = 0.0;
    }
    // end if alphaAverage
  }
  // end if ibs coupling
  coeffMulT = sigs * 2.0 * sqrt(pi) / (bunchParam["nMacro"] * dtsamp2 * clight);

  std::vector<double> ibscoeffs;
  ibscoeffs.push_back(coeffx);
  ibscoeffs.push_back(coeffy);
  ibscoeffs.push_back(coeffs);
  ibscoeffs.push_back(coeffMulT);

  ibsCoeff.push_back(ibscoeffs);
}

void Bunch::printBunchParameters() {
  std::printf("%-30s\n", "Bunch Parameters");
  std::printf("%-30s\n", "================");
  for (auto const &pair : bunchParam) {
    std::printf("%-30s %16.8e\n", pair.first.c_str(), pair.second);
  }
  std::printf("\n");
}

void Bunch::printTwissHeader() {
  std::printf("%-30s\n", "Twiss Header");
  std::printf("%-30s\n", "============");
  for (auto const &pair : twheader) {
    std::printf("%-30s %16.8e\n", pair.first.c_str(), pair.second);
  }
  std::printf("\n");
}

void Bunch::printLongParam() {
  std::printf("%-30s\n", "Longitudinal Parameters");
  std::printf("%-30s\n", "=======================");
  for (auto const &pair : longitudinalParameters) {
    std::printf("%-30s %16.8e\n", pair.first.c_str(), pair.second);
  }
  std::printf("\n");
}
void Bunch::printRadiationParameters() {
  std::printf("%-30s\n", "Radiation Parameters");
  std::printf("%-30s\n", "===================");
  for (auto const &pair : Bunch::radiationParameters) {
    std::printf("%-30s %16.8e\n", pair.first.c_str(), pair.second);
  }
  std::printf("\n");
}

void Bunch::printDistribution() {
  std::for_each(distribution.begin(), distribution.end(),
                [](std::vector<double> &particle) {
                  std::printf("%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
                              particle[0], particle[1], particle[2],
                              particle[3], particle[4], particle[5]);
                });

  std::printf("\n");
};

void Bunch::printEmittance() {
  std::for_each(emittances.begin(), emittances.end(),
                [](std::vector<double> &particle) {
                  std::printf("%12.8e %12.8e %12.8e %12.8e\n", particle[0],
                              particle[1], particle[2], particle[3]);
                });

  std::printf("\n");
};

void Bunch::printIBSGrowthRates() {
  std::for_each(ibsGrowthRates.begin(), ibsGrowthRates.end(),
                [](std::vector<double> &particle) {
                  std::printf("%12.8e %12.8e %12.8e\n", particle[0],
                              particle[1], particle[2]);
                });

  std::printf("\n");
};

#define MAX_DATE 12

std::string Bunch::get_date() {

  /* function returning current date as string for creating timestamped output
   * files */

  time_t now;

  // struct tm * timeinfo;
  char the_date[MAX_DATE];
  // char buffer [80];

  // time_t rawtime;
  // timeinfo = localtime (&rawtime);
  // strftime (buffer,80,"Now it's %I:%M%p.",timeinfo);

  // std::stringstream ss;
  // ss << fn << "_" << "%d_%m_%Y" << "_" << buffer << ".dat";
  // std::string s = ss.str();

  the_date[0] = '\0';

  now = time(NULL);

  if (now != -1) {
    strftime(the_date, MAX_DATE, "%d_%m_%Y", gmtime(&now));
  }

  return std::string(the_date);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (int i = 0; i < v.size(); i++) {
    os << std::scientific << std::setw(15) << v[i];
  }
  // os << std::endl;
  // std::copy(v.begin(), v.end(), ostream_iterator<T>(os, "\t"));
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os,
                         const std::vector<std::vector<T>> &v) {
  using namespace std;

  // NOTE: for some reason std::copy doesn't work here, so I use manual loop
  // copy(v.begin(), v.end(), ostream_iterator<std::vector<T>>(os, "\n"));

  for (size_t i = 0; i < v.size(); ++i)
    os << std::setw(6) << i << v[i] << "\n";
  return os;
}

void Bunch::writeDistribution(int turn) {

  std::stringstream ss;
  ss << "Distribution_bucket_" << bunchParam["bucket"] << "_"
     << get_date()
     // << "_turn_" << std::setw(10) << std::setfill('0') << turn
     << ".dat";

  std::string s;
  s = ss.str();

  if (turn == 0) {
    if (remove(s.c_str()) != 0)
      perror("Error deleting file");
    else
      puts("File successfully deleted");
  }
  std::ofstream ofile;
  ofile.open(s.c_str(), std::ios_base::app);
  ofile << distribution;
  // std::copy(distribution.begin(), distribution.end(),
  //        std::ostream_iterator<std::vector<double>>(ofile));

  // std::cout << "writing vector size = " << distribution.size() << std::endl;
  ofile.close();
}

void Bunch::writeEmittances() {

  std::stringstream ss;
  ss << "CTE_Emittances_" << bunchParam["bucket"] << "_" << get_date()
     << ".dat";

  std::string s;
  s = ss.str();

  if (remove(s.c_str()) != 0)
    perror("Error deleting file");
  else
    puts("File successfully deleted");

  std::ofstream ofile(s.c_str());
  ofile << emittances;

  ofile.close();
}