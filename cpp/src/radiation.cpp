#include "../include/cte_bits/random.hpp"
#include <algorithm>
#include <complex>
#include <ibs>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>

namespace cte_radiation {
/*
 * =============================================================================
 * =============================================================================
 * GET RADIATION DAMPING RELEVANT VALUES.
 * =============================================================================
 * AUTHORS:
 *  - TOM MERTENS
 *
 * REFS:
 * =============================================================================
 * Arguments:
 * ----------
 *      - std::map<std::string, double> &twheader
 *          MADX twiss summary (read with IBSLib GetTwissHeader)
 * Returns:
 * --------
 *      - std::map<std::string, double>
 *          output["taux"] = 1.0 / alphax;
 *          output["tauy"] = 1.0 / alphay;
 *          output["taus"] = 1.0 / alphas;
 *          output["exinf"] = exinf;
 *          output["eyinf"] = eyinf;
 *          output["sigeoe2"] = sigE0E2;
 *          output["sigsinf"] = sigs;
 *          output["jx"] = jx;
 *          output["jy"] = jy;
 * =============================================================================
 * =============================================================================
 */
std::map<std::string, double>
radiationEquilib(std::map<std::string, double> &twheader,
                 std::map<std::string, double> &longparam) {
  const double electron_volt_joule_relationship = 1.602176634e-19;
  const double hbar = 1.0545718176461565e-34;
  std::map<std::string, double> output;

  // base quantities
  double gamma = twheader["GAMMA"];
  double gammatr = twheader["GAMMATR"];
  double p0 = twheader["PC"] * 1.0e9;
  double len = twheader["LENGTH"] * 1.0;
  double restE = twheader["MASS"] * 1.0e9;
  double charge = twheader["CHARGE"] * 1.0;
  double aatom = twheader["aatom"];
  double qs = longparam["qs"];
  double omega = twheader["omega"];
  double q1 = twheader["Q1"];

  // use madx rad integrals
  double I1 = twheader["SYNCH_1"];
  double I2 = twheader["SYNCH_2"];
  double I3 = twheader["SYNCH_3"];
  double I4x = twheader["SYNCH_4"];
  double I5x = twheader["SYNCH_5"];

  double I4y = 0.0;
  double I5y = 0.0;

  // derived quantities
  double pradius = ParticleRadius(charge, aatom);
  double CalphaEC =
      pradius * clight / (3.0 * restE * restE * restE) * (p0 * p0 * p0 / len);

  // transverse partition numbers
  double jx = 1.0 - I4x / I2;
  double jy = 1.0 - I4y / I2;
  double alphax = 2.0 * CalphaEC * I2 * jx;
  double alphay = 2.0 * CalphaEC * I2 * jy;
  double alphas = 2.0 * CalphaEC * I2 * (jx + jy);

  // mc**2 expressed in Joule to match units of cq
  double mass = restE * electron_volt_joule_relationship;
  double cq = 55.0 / (32.0 * sqrt(3.0)) * (hbar * clight) / mass;

  double sigE0E2 = cq * gamma * gamma * I3 / (2.0 * I2 + I4x + I4y);
  // ! = deltaE/E_0 see wiedemann p. 302,
  // and Wolski: E/(p0*c) - 1/beta0 = (E - E0)/(p0*c) = \Delta E/E0*beta0 with
  // E0 = p0*c/beta0 therefore:
  double betar = BetaRelativisticFromGamma(gamma);
  double dpop = dee_to_dpp(sqrt(sigE0E2), betar);
  double sigs = dpop * len * eta(gamma, gammatr) / (2 * pi * qs);
  double exinf = cq * gamma * gamma * I5x / (jx * I2);
  double eyinf = cq * gamma * gamma * I5y / (jy * I2);

  double betaAvg = len / (q1 * 2.0 * pi);

  eyinf = (eyinf == 0.0) ? cq * betaAvg * I3 / (2.0 * jy * I2) : eyinf;

  output["taux"] = 1.0 / alphax;
  output["tauy"] = 1.0 / alphay;
  output["taus"] = 1.0 / alphas;
  output["exinf"] = exinf;
  output["eyinf"] = eyinf;
  output["sigeoe2"] = sigE0E2;
  output["sigsinf"] = sigs;
  output["jx"] = jx;
  output["jy"] = jy;

  return output;
}

double RadiationLossesPerTurn(std::map<std::string, double> &twiss) {
  double gamma = twiss["GAMMA"];
  double p0 = twiss["PC"];
  double len = twiss["LENGTH"];
  double mass = twiss["MASS"];
  double charge = twiss["CHARGE"];
  double aatom = twiss["aatom"];
  double I2 = twiss["SYNCH_2"];

  double particle_radius = ParticleRadius(charge, aatom);
  double cgamma = (4.0 * pi / 3.0) * (particle_radius / (mass * mass * mass));
  double trev = twiss["trev"];

  return (clight * cgamma) / (2.0 * pi * len) * p0 * p0 * p0 * p0 * I2 * 1.0e9 *
         trev;
}

void CalcRadDecayExcitation(std::map<std::string, double> &twheader,
                            std::map<std::string, double> &radparam) {
  double gamma = twheader["GAMMA"];
  double trev = twheader["trev"];
  double timeratio = twheader["timeratio"];
  double taus = radparam["taus"];
  double taux = radparam["taux"];
  double tauy = radparam["tauy"];
  double sigsinf = radparam["sigsinf"];
  double tt = trev * timeratio;
  double sigx = sqrt(radparam["exinf"] * twheader["betxavg"]);
  double sigy = sqrt(radparam["eyinf"] * twheader["betyavg"]);

  // timeratio is real machine turns over per simulation turn
  radparam["coeffdecaylong"] = exp(-tt / taus);

  // excitation uses a uniform distibution on [-1:1]
  // sqrt(3) * sigma => +/-3 sigma**2
  // see also lecture 2 Wolski on linear dynamics and radiation damping
  // not sure about the sqrt3
  radparam["coeffexcitelong"] = (radparam["sigeoe2"] * gamma) / sqrt(3) *
                                sqrt(3.) * sqrt(1.0 * tt / taus);

  // the damping time is for EMITTANCE, therefore need to multiply by 2
  radparam["coeffdecayx"] = exp(-(tt / (2 * taux)));
  radparam["coeffdecayy"] = exp(-(tt / (2 * tauy)));

  // exact:
  // coeffgrow= sigperp*sqrt(3.)*sqrt(1-coeffdecay**2)
  // squared because sigma and not emit
  radparam["coeffgrowx"] =
      sigx * sqrt(3.) * sqrt(1 - pow(radparam["coeffdecayx"], 2));
  radparam["coeffgrowy"] =
      sigy * sqrt(3.) * sqrt(1 - pow(radparam["coeffdecayy"], 2));
};

} // namespace cte_radiation