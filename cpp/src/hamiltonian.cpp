#include <algorithm>
#include <ibs>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace cte_hamiltonian {
/*
================================================================================
================================================================================
METHOD TO CALCULATE TIME COEFFICIENT.
================================================================================
AUTHORS:
    -  TOM MERTENS

HISTORY:
    - initial version 02/08/2021

REF:
    - based on original CTE code (CERN) - but updated for multi-RF systems
    - Lee (Third Ed.) page 233 eq. 3.13 (kinetic part of Hamiltonian)
================================================================================
Arguments:
----------
    - std::map<std::string, double> &twissheaderL
        twissheader map updated with longitudinal parameters
        ( ref: updateTwissHeaderLong )
    - double baseHarmonicNumber
        the base harmonic number of the main RF freq (defining the bucket
        number)

Returns:
--------
    - double tcoeff

================================================================================
================================================================================
*/
double tcoeff(std::map<std::string, double> &twissheaderL,
              double baseHarmonicNumber) {
  return twissheaderL["omega"] * twissheaderL["eta"] * baseHarmonicNumber;
}

/*
================================================================================
================================================================================
METHOD TO CALCULATE MOMENTUM COEFFICIENT.
================================================================================
AUTHORS:
    -  TOM MERTENS

HISTORY:
    - initial version 02/08/2021

REF:
    - based on original CTE code (CERN) - but updated for multi-RF systems
    - PCOEFF is the coefficient in front of the goniometric functions
      (potential) in the Hamiltonian
    - Lee (Third Ed.) page 233 eq. 3.13

================================================================================
Arguments:
----------
    - std::map<std::string, double> &twissheaderL
        twissheader map updated with longitudinal parameters
        ( ref: updateTwissHeaderLong )

Returns:
--------
    - double pcoeff

================================================================================
================================================================================
*/
double pcoeff(std::map<std::string, double> &twissheaderL, double voltage) {
  // factor 1.0e9 -> pc is in GeV
  return twissheaderL["omega"] * voltage * twissheaderL["CHARGE"] /
         (2.0 * pi * twissheaderL["PC"] * 1.0e9 * twissheaderL["betar"]);
}
/*
================================================================================
================================================================================
METHOD TO CALCULATE (APPROX) LONGITUDINAL HAMILTONIAN IN FUNCTION OF GIVEN t
(for synchronuous particle t=0).
================================================================================
AUTHORS:
    -  TOM MERTENS

HISTORY:
    - initial version 02/08/2021

REF:
    - based on original CTE code (CERN) - but updated for multi-RF systems
    - Lee (Third Ed.) page 233 eq. 3.13


HAMILTONIAN:

H = 0.5 * tcoeff * delta**2 +
SUM_i(pcoeff[v_i,omega0,p0,betarelativistic](cos(h_i omega0 t)- cos(h_i/h_0
phis) + (h_i omega0 t- h_i/h_0 phis) sin(h_i/h_0 phis)))

================================================================================
Arguments:
----------
    - std::map<std::string, double> &twissheaderL
        twissheader map updated with longitudinal parameters
        ( ref: updateTwissHeaderLong )
    - std::vector<double> &harmonicNumbers
        RF harmonic numbers - Assuming main harmonic number is first in the list
    - std::vector<double> &rfvoltages
        RF voltages
    - double tcoeff
        see tcoeff
    - double t
        time point to calculate Hamiltonian (converted to phase internally)

Returns:
--------
    - double Hamiltonian value

================================================================================
================================================================================
*/
double Hamiltonian(std::map<std::string, double> &twheader,
                   std::map<std::string, double> &longparam,
                   std::vector<double> &h, std::vector<double> &v,
                   double tcoeff, double t, double delta) {
  double kinetic, potential;

  // kinetic contribution
  // We assume initial bunch length is given
  kinetic = 0.5 * tcoeff * delta * delta;

  std::vector<double> pcoeffs, hRatios, hRatiosInv, phases;

  // calculate coefficients for the determining the potential
  for (int i = 0; i < h.size(); i++) {
    pcoeffs.push_back(pcoeff(twheader, v[i]));
    phases.push_back(h[i] * twheader["omega"] * t);
    hRatios.push_back(h[0] / h[i]);
    hRatiosInv.push_back(h[i] / h[0]);
  }

  // calc the potential
  potential =
      pcoeffs[0] * (cos(phases[0]) - cos(longparam["phis"]) +
                    (phases[0] - longparam["phis"]) * sin(longparam["phis"]));

  for (int i = 1; i < h.size(); i++) {
    potential += pcoeffs[i] * hRatios[i] *
                 (cos(phases[i]) - cos(hRatiosInv[i] * longparam["phis"]) +
                  (phases[i] - hRatiosInv[i] * longparam["phis"]) *
                      sin(hRatiosInv[i] * longparam["phis"]));
  }

  return kinetic + potential;
}
} // namespace cte_hamiltonian