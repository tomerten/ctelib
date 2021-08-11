#include <cte>
#include <ibs>
#include <map>
#include <numeric>
#include <string>
#include <vector>
int main() {
  string twissfilename = "../src/b2_design_lattice_1996.twiss";
  std::map<std::string, double> twheader;
  twheader = GetTwissHeader(twissfilename);
  std::map<std::string, std::vector<double>> tw =
      GetTwissTableAsMap(twissfilename);

  // rf settings
  std::vector<double> h, v;
  h.push_back(400.0);
  v.push_back(-1.5e6);

  std::map<std::string, double> bp;
  bp["bucket"] = 0;
  bp["atomNumber"] = emass / pmass;
  bp["charge"] = -1;
  bp["nMacro"] = 100;
  bp["nReal"] = 1.e10;
  bp["sigs"] = 0.005;
  bp["seed"] = 123456;
  bp["ex"] = 5.5e-9;
  bp["ey"] = 0.05 * bp["ex"];
  bp["timeRatio"] = 1.0;

  Bunch testbunch(twheader, tw, bp, h, v);
  testbunch.printBunchParameters();
  testbunch.printTwissHeader();
  testbunch.printLongParam();
  testbunch.printRadiationParameters();
  testbunch.printDistribution();
  // testbunch.writeDistribution(0);

  double K2L = std::accumulate(tw["K2L"].begin(), tw["K2L"].begin(), 0.0);
  double K2SL = std::accumulate(tw["K2SL"].begin(), tw["K2SL"].begin(), 0.0);

  for (int turn = 0; turn < 10; turn++) {
    PHYSICS::updateBeta(testbunch.distribution, testbunch.twheader, 0.0, K2L,
                        K2SL);

    PHYSICS::updateRad(testbunch.distribution, testbunch.radiationParameters,
                       bp["seed"]);
    PHYSICS::updateRF(testbunch.distribution, 1.0, testbunch.twheader,
                      testbunch.longitudinalParameters, h, v);

    if (int(turn % 10) == 0) {
      testbunch.writeDistribution(turn);
      testbunch.getEmittance();
    }
  }
  testbunch.printEmittance();
  testbunch.writeEmittances();
  testbunch.getIBSGrowthRates(13);
  testbunch.printIBSGrowthRates();

  ParticlesTimesToHistogram(testbunch.distribution, 100,
                            testbunch.longitudinalParameters["tauhat"] / 2.,
                            testbunch.longitudinalParameters["phisNext"] /
                                (h[0] * testbunch.twheader["omega"]));
  return 0;
}