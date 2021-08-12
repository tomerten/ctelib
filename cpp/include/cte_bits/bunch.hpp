#include <map>
#include <string>
#include <vector>

// TODO make this a map so one can give a dict in python
struct bunchParameters {
  // bucket number
  int bucket;
  // atomic mass number
  double atomNumber;
  // particle charge
  int charge;
  // number of macro particle to use
  int nMacro;
  // real number of particles
  double nReal;
  // bunch length
  double sigs;
  // random seed
  int seed;
};

class Bunch {
public:
  std::map<std::string, double> bunchParam;
  std::map<std::string, std::vector<double>> tw;
  std::map<std::string, double> twheader;
  std::vector<std::vector<double>> distribution;
  std::vector<std::vector<double>> emittances;
  std::vector<std::vector<double>> ibsGrowthRates;
  std::vector<std::vector<double>> ibsCoeff;
  std::map<std::string, std::vector<double>> bunchParametersLocal;
  std::map<std::string, double> radiationParameters;
  std::map<std::string, double> longitudinalParameters;
  // constructor
  Bunch(std::map<std::string, double> &,
        std::map<std::string, std::vector<double>> &,
        std::map<std::string, double> &, std::vector<double> &,
        std::vector<double> &);

  // helper functions
  void resetDebunchLosses();
  std::string get_date();
  void getEmittance();
  void getIBSGrowthRates(int);
  void getIBSCoefficients();

  // set functions
  void setBasic();
  void setDistribution(std::vector<double> &h, std::vector<double> &v);
  void setLongitudinalParameters(std::vector<double> &, std::vector<double> &);
  void setRadiationParameters();

  // print to screen methods
  void printBunchParameters();
  void printTwissHeader();
  void printLongParam();
  void printRadiationParameters();
  void printDistribution();
  void printEmittance();
  void printIBSGrowthRates();
  void printHistogramTime();
  void printSqrtHistorgram();

  // write to file methods
  void writeDistribution(int);
  void writeEmittances();
  void writeLocalParameters();
};