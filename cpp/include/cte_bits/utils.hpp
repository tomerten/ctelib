#include <vector>

struct ParticleTimesToInteger {
  /*
  0 -> tauahat = half the phase acceptance
  1 -> nbins = numbers of bins to use in binning
  2 -> t synchronous to shift the center of the distribution to zero for the
  binning
  */
  int operator()(std::vector<double> &data, std::vector<double> &params);
};
void denseHistogram(std::vector<int> data, std::vector<int> &histogram,
                    int nbins, double tauhat);
std::vector<int>
ParticlesTimesToHistogram(std::vector<std::vector<double>> &data, int nbins,
                          double tauhat, double ts);
std::vector<double> getColumnMeans(std::vector<std::vector<double>> &dist);
std::vector<double> vectorMultiply(std::vector<double> x,
                                   std::vector<double> y);

std::vector<double> vectorAdd(std::vector<double> x, std::vector<double> y);

std::vector<double> vectorSub(std::vector<double> x, std::vector<double> y);
std::vector<double> CalcRMS(const std::vector<std::vector<double>> dist);