#include <chrono>
#include <cte>
#include <ibs>
#include <map>
#include <numeric>
#include <string>
#include <vector>
int main() {
  auto t1 = std::chrono::high_resolution_clock::now();
  string twissfilename = "../src/b2_design_lattice_1996.twiss";
  std::map<std::string, double> twheader;
  twheader = GetTwissHeader(twissfilename);
  std::map<std::string, std::vector<double>> tw =
      GetTwissTableAsMap(twissfilename);

  // rf settings
  std::vector<double> h, v;
  h.push_back(400.0);
  v.push_back(-1.5e6);

  std::map<std::string, double> bp, in;
  bp["bucket"] = 0;
  bp["atomNumber"] = emass / pmass;
  bp["charge"] = -1;
  bp["nMacro"] = 1000;
  bp["nReal"] = 1.e10;
  bp["nbins"] = 100;
  bp["sigs"] = 0.005;
  bp["seed"] = 123456;
  bp["ex"] = 5.5e-9;
  bp["ey"] = 0.05 * bp["ex"];
  bp["timeRatio"] = 100.0;
  bp["fracibstot"] = 1.0;
  bp["ibsCoupling"] = 0.0;
  bp["model"] = 1;
  in = readInput("testinput.in");

  Bunch testbunch(twheader, tw, in, h, v);
  testbunch.printBunchParameters();
  testbunch.printTwissHeader();
  testbunch.printLongParam();
  testbunch.printRadiationParameters();
  // testbunch.printDistribution();
  // testbunch.writeDistribution(0);

  double K2L = std::accumulate(tw["K2L"].begin(), tw["K2L"].begin(), 0.0);
  double K2SL = std::accumulate(tw["K2SL"].begin(), tw["K2SL"].begin(), 0.0);
  testbunch.getEmittance();

  int barWidth = 70;
  int nturns = (int)in["nturns"];
  int nwrite = (int)in["nwrite"];

  std::printf("nturns %i", nturns);

  for (int turn = 0; turn < nturns; turn++) {

    // progress bar
    std::cout << "[";
    int progress = (double)turn / nturns * barWidth;
    for (int j = 0; j < barWidth; ++j) {
      if (j < progress)
        std::cout << "=";
      else if (j == progress)
        std::cout << ">";
      else
        std::cout << " ";
    }
    std::cout << "]" << int((double)turn / nturns * 100) << " %\r";
    std::cout.flush();

    PHYSICS::updateBeta(testbunch.distribution, testbunch.twheader, 0.0, K2L,
                        K2SL);

    PHYSICS::updateRad(testbunch.distribution, testbunch.radiationParameters,
                       bp["seed"]);
    PHYSICS::updateRF(testbunch.distribution, 1.0, testbunch.twheader,
                      testbunch.longitudinalParameters, h, v,
                      testbunch.debunchLosses);

    testbunch.updateIBS(h);
    PHYSICS::updateIBS(testbunch.distribution, testbunch.longitudinalParameters,
                       testbunch.bunchParam, testbunch.twheader, h,
                       testbunch.ibsCoeff.back(), testbunch.sqrthistogram);
    testbunch.getEmittance();

    if (int(turn % nwrite) == 0) {
      testbunch.writeDistribution(turn);
    }
  }
  std::cout << std::endl;
  testbunch.printDebunchLosses();

  // testbunch.printEmittance();
  testbunch.writeEmittances();
  // testbunch.printIBSGrow // end progressbar
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count() / 1000.0 << "s";
  return 0;
}