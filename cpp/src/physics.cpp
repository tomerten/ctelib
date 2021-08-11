
#include "../include/cte_bits/hamiltonian.hpp"
#include "../include/cte_bits/random.hpp"
#include "../include/cte_bits/synch.hpp"
#include <algorithm>
#include <complex>
#include <ibs>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>

namespace PHYSICS {
std::vector<double> vectorMultiply(std::vector<double> x,
                                   std::vector<double> y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(),
                 std::multiplies<double>());
  return x;
}

std::vector<double> vectorAdd(std::vector<double> x, std::vector<double> y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<double>());
  return x;
}

// apply radiation damping
void updateRad(std::vector<std::vector<double>> &dist,
               std::map<std::string, double> &radparam, int seed) {

  // create decay vector 6D
  std::vector<double> D;
  D.push_back(radparam["coeffdecayx"]);
  D.push_back(radparam["coeffdecayx"]);
  D.push_back(radparam["coeffdecayy"]);
  D.push_back(radparam["coeffdecayy"]);
  D.push_back(1.0);
  D.push_back(radparam["coeffdecaylong"]);

  // create excitation vector 6D
  std::vector<double> E;
  E.push_back(radparam["coeffgrowx"]);
  E.push_back(radparam["coeffgrowx"]);
  E.push_back(radparam["coeffgrowy"]);
  E.push_back(radparam["coeffgrowy"]);
  E.push_back(0.0);
  E.push_back(radparam["coeffexcitelong"]);

  // generate vector of const vectors
  // for decay, excitation and random coefficients
  int n = dist.size();
  std::vector<std::vector<double>> decay(n);
  std::vector<std::vector<double>> excitation(n);
  std::vector<std::vector<double>> R(n);

  std::fill(decay.begin(), decay.end(), D);
  std::fill(excitation.begin(), excitation.end(), E);

  // generate raomdom vectors
  std::for_each(R.begin(), R.end(), [&seed](std::vector<double> &particle) {
    for (int i = 0; i < 6; i++) {
      particle.push_back(2.0 * cte_random::ran3(&seed) - 1.0);
    }
  });

  // update the distribution
  // dist * decay + excitation * random

  // dist * decay
  std::transform(dist.begin(), dist.end(), decay.begin(), dist.begin(),
                 vectorMultiply);

  // excitation * random
  std::transform(excitation.begin(), excitation.end(), R.begin(),
                 excitation.begin(), vectorMultiply);

  // dist * decay + excitation * random
  std::transform(dist.begin(), dist.end(), excitation.begin(), dist.begin(),
                 vectorAdd);

  // print out for debug
  /*
  std::for_each(dist.begin(), dist.end(), [](std::vector<double> &particle) {
    std::printf("%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", particle[0],
                particle[1], particle[2], particle[3], particle[4],
                particle[5]);
  });

  std::printf("\n");
  */
}

void updateRF(std::vector<std::vector<double>> &dist, double fmix,
              std::map<std::string, double> &tw,
              std::map<std::string, double> &longparam, std::vector<double> &hs,
              std::vector<double> &vs) {

  std::for_each(
      dist.begin(), dist.end(),
      [&tw, &hs, &fmix, &vs, &longparam](std::vector<double> &particle) {
        double charge = tw["CHARGE"];
        double phis = longparam["phis"];
        double gamma = tw["GAMMA"];
        double p0 = tw["PC"] * 1.0e9;
        double betar = BetaRelativisticFromGamma(gamma);

        double omega = tw["omega"];
        double h0 = hs[0];
        double tc = cte_hamiltonian::tcoeff(tw, h0);

        // Lee Third edition page 233 eqs 3.6 3.16 3.13
        // the phase  = h0 omega0 t
        double phi = h0 * omega * particle[4];

        phi += fmix * 2.0 * pi * h0 * tw["eta"] * particle[5];
        particle[4] = phi / (h0 * omega);

        double voltdiff = CTESYNCH::VoltageRf(phi, vs, hs);
        voltdiff -= CTESYNCH::VoltageRf(phis, vs, hs);

        particle[5] += fmix * voltdiff * charge / (betar * betar * p0);
      });

  // print out for debug
  /*
  std::for_each(dist.begin(), dist.end(), [](std::vector<double> &particle) {
    std::printf("%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", particle[0],
                particle[1], particle[2], particle[3], particle[4],
                particle[5]);
  });

  std::printf("\n");
  */
}

void updateBeta(std::vector<std::vector<double>> &dist,
                std::map<std::string, double> tw, double coupling, double K2L,
                double K2SL) {

  double qx = tw["Q1"];
  double qy = tw["Q2"];
  double ksix = tw["DQ1"];
  double ksiy = tw["DQ2"];

  double psix = 2.0 * pi * qx;
  double psiy = 2.0 * pi * qy;
  std::for_each(dist.begin(), dist.end(),
                [&tw, &psix, &ksix, &psiy, &ksiy, &coupling, &K2L,
                 &K2SL](std::vector<double> &particle) {
                  // rotation in x-px plane
                  double psi1 = psix + particle[5] * ksix;
                  double a11 = cos(psi1);
                  double a12 = sin(psi1);

                  double temp = particle[0];
                  particle[0] = particle[0] * a11 + particle[1] * a12;
                  particle[1] = particle[1] * a11 - temp * a12;

                  // rotation in y-py plane
                  double psi2 = psiy + particle[5] * ksiy;
                  a11 = cos(psi2);
                  a12 = sin(psi2);

                  temp = particle[2];
                  particle[2] = particle[2] * a11 + particle[3] * a12;
                  particle[3] = particle[3] * a11 - temp * a12;

                  // now have dqmin part - coupling between x and y
                  particle[1] += coupling * particle[2];
                  particle[3] += coupling * particle[0];

                  // thin sextupole kick
                  particle[1] +=
                      0.5 * K2L * (pow(particle[0], 2) - pow(particle[2], 2)) -
                      K2SL * (particle[0] * particle[2]);
                  particle[3] +=
                      0.5 * K2SL * (pow(particle[0], 2) - pow(particle[2], 2)) +
                      K2L * (particle[0] * particle[2]);
                });
}

} // namespace PHYSICS