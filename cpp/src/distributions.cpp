#include "../include/cte_bits/hamiltonian.hpp"
#include "../include/cte_bits/random.hpp"
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <cstdlib>
#include <ibs>
#include <iostream>
#include <map>
#include <math.h>
#include <string>

namespace cte_distributions {
std::vector<double> BiGaussian4D(const double betax, const double ex,
                                 const double betay, const double ey,
                                 int seed) {
  std::vector<double> out;

  // std::printf("%-30s %12.6e \n", "betx", betax);
  // std::printf("%-30s %12.6e \n", "bety", betay);
  // std::printf("%-30s %12.6e \n", "ex", ex);
  // std::printf("%-30s %12.6e \n", "ey", ey);
  static double ampx, ampy, amp, r1, r2, facc;
  static double x, px, y, py;

  // 1 sigma rms beam sizes using average ring betas
  ampx = sqrt(betax * ex);
  ampy = sqrt(betay * ey);

  // generate bi-gaussian distribution in the x-px phase-space
  do {
    r1 = 2 * cte_random::ran3(&seed) - 1;
    r2 = 2 * cte_random::ran3(&seed) - 1;
    amp = r1 * r1 + r2 * r2;
  } while (amp >= 1);

  facc =
      sqrt(-2 * log(amp) /
           amp); // transforming [-1,1] uniform to gaussian - inverse transform

  x = ampx * r1 * facc;  // scaling the gaussian
  px = ampx * r2 * facc; // scaling the gaussian

  // generate bi-gaussian distribution in the y-py phase-space
  do {
    r1 = 2 * cte_random::ran3(&seed) - 1;
    r2 = 2 * cte_random::ran3(&seed) - 1;
    amp = r1 * r1 + r2 * r2;
  } while (amp >= 1);

  facc =
      sqrt(-2 * log(amp) /
           amp); // transforming [-1,1] uniform to gaussian - inverse transform

  y = ampy * r1 * facc;  // scaling the gaussian
  py = ampy * r2 * facc; // scaling the gaussian

  // std::printf("%12.6e %12.6e %12.6e %12.6e\n", x, px, y, py);
  out.push_back(x);
  out.push_back(px);
  out.push_back(y);
  out.push_back(py);

  return out;
}

std::vector<double>
BiGaussian6DLongMatched(double betax, double ex, double betay, double ey,
                        std::vector<double> &h, std::vector<double> &v,
                        std::map<std::string, double> &twheader,
                        std::map<std::string, double> &longparam, int seed) {

  double h0 = h[0];
  double tauhat = twheader["tauhat"];
  double omega = twheader["omega"];
  double ampt = longparam["sigs"] / clight;
  // Max value Hamiltonian that is stable
  // is, with the sign convention used, left of the ham contour
  // at 180-phis (The Ham rises lin to the right.)
  double hom = (h0 * omega);
  double ts = longparam["phis"] / hom;
  int npi = int(longparam["phis"] / (2.0 * pi));
  double tperiod = 2.0 * pi / hom;
  double ts2 = longparam["phisNext"] / hom + double(npi) * tperiod;
  double delta = longparam["sige"];

  std::vector<double> out;
  out = BiGaussian4D(betax, ex, betay, ey, seed);

  // adding two zeros
  out.push_back(0.0);
  out.push_back(0.0);

  double r1, r2, amp, facc, tc, pc, ham, hammin;
  tc = (omega * twheader["eta"] * h0);
  pc = (omega * twheader["CHARGE"]) /
       (2.0 * pi * twheader["PC"] * 1.0e9 * twheader["betar"]);

  // std::printf("%-20s %16.8e\n", "tc", tc);
  // std::printf("%-20s %16.8e\n", "pc", pc);
  // max Hamiltonian
  double hammax =
      cte_hamiltonian::Hamiltonian(twheader, longparam, h, v, tc, ts, 0.0);
  /*
  std::printf("%-20s %16.8e\n", "hammax", hammax);
  std::printf("%-20s %16.8e\n", "ts", ts);
  std::printf("%-20s %16.8e\n", "ts2", ts2);
  */
  // select valid t values
  do {
    // looper++;
    // std::printf("%i\n", looper);
    r1 = 2 * cte_random::ran3(&seed) - 1;
    r2 = 2 * cte_random::ran3(&seed) - 1;
    amp = r1 * r1 + r2 * r2;
    if (amp >= 1)
      continue;

    facc = sqrt(-2 * log(amp) / amp);
    // std::printf("%-20s %16.8e\n", "out4", out[4]);
    out[4] = ts2 + ampt * r1 * facc;
    /*
        std::printf("%-20s %16.8e\n", "ts", ts);
        std::printf("%-20s %16.8e\n", "ampt", ampt);
        std::printf("%-20s %16.8e\n", "r1", r1);
        std::printf("%-20s %16.8e\n", "facc", facc);
        std::printf("%-20s %16.8e\n", "out4", out[4]);
        std::printf("%-20s %16.8e\n", "out4-ts", out[4] - ts);
        std::printf("%-20s %16.8e\n", "tauhat", tauhat);
    */
    if (abs(out[4] - ts2) >= abs(ts - ts2))
      continue;

    // min Hamiltonian
    hammin = cte_hamiltonian::Hamiltonian(twheader, longparam, h, v, tc, out[4],
                                          0.0);
    // std::printf("%-20s %16.8e\n", "hammin", hammin);
    // std::printf("%-20s %16.8e\n", "hammax", hammax);

  } while ((hammin > hammax) || (abs(out[4] - ts2) >= abs(ts - ts2)));

  // select matched deltas
  do {
    // looper++;
    // std::printf("%i\n", looper);
    r1 = 2 * cte_random::ran3(&seed) - 1;
    r2 = 2 * cte_random::ran3(&seed) - 1;
    amp = r1 * r1 + r2 * r2;

    if (amp >= 1)
      continue;

    facc = sqrt(-2 * log(amp) / amp);
    out[5] = longparam["sige"] * r2 * facc;
    // std::printf("%-20s %16.8e\n", "delta", out[5]);

    ham = cte_hamiltonian::Hamiltonian(twheader, longparam, h, v, tc, out[4],
                                       out[5]);
    // std::printf("%-20s %16.8e\n", "hammin", hammin);
    // std::printf("%-20s %16.8e\n", "ham", ham);
    // std::printf("%-20s %16.8e\n", "hammax", hammax);
  } while ((ham < hammin) || (ham > hammax));

  return out;
}

std::vector<std::vector<double>> GenerateDistributionMatched(
    int nMacro, double betax, double ex, double betay, double ey,
    std::vector<double> &h, std::vector<double> &v,
    std::map<std::string, double> &twheader,
    std::map<std::string, double> &longparam, int seed) {
  std::vector<std::vector<double>> out;

  for (int i = 0; i < nMacro; i++) {
    out.push_back(BiGaussian6DLongMatched(betax, ex, betay, ey, h, v, twheader,
                                          longparam, seed));
  }
  return out;
}

} // namespace cte_distributions