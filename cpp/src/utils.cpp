#include <algorithm>
//#include <boost/iterator/counting_iterator.hpp>
#include <ibs>
#include <iterator>
#include <math.h>
#include <numeric>
#include <vector>
// gets the particle times and uses the phase acceptance and number of bins
// to generate an integer which is than binned/histogram by
// ParticleTimesToHistogram
struct ParticleTimesToInteger {
  /*
  0 -> tauahat = half the phase acceptance
  1 -> nbins = numbers of bins to use in binning
  2 -> t synchronous to shift the center of the distribution to zero for the
  binning
  */
  int operator()(std::vector<double> &data, std::vector<double> &params) const {
    double out;
    double dtsamp2 = 2 * params[0] / params[1];
    out = (data[4] - params[2] + params[0]) / dtsamp2;
    out = (int)(out + 0.5f);
    return out;
  }
};

void denseHistogram(std::vector<int> data, std::vector<int> &histogram,
                    int nbins, double tauhat) {
  std::vector<int> input = data;
  std::sort(input.begin(), input.end());

  int numBins = nbins + 1;
  histogram.resize(numBins);
  for (int i = 0; i < numBins; i++) {
    auto upper = std::upper_bound(input.begin(), input.end(), i);
    histogram[i] = std::distance(input.begin(), upper);
  }
  /*
    for (std::vector<int>::iterator i = histogram.begin(); i != histogram.end();
         i++) {
      std::printf("%12i\n", *i);
    }
  */
  std::adjacent_difference(histogram.begin(), histogram.end(),
                           histogram.begin());
}

std::vector<int>
ParticlesTimesToHistogram(std::vector<std::vector<double>> &data, int nbins,
                          double tauhat, double ts) {

  int n = data.size();
  std::vector<int> timecomponent(n);
  std::vector<int> histogram;
  std::vector<double> param;

  param.push_back(tauhat);
  param.push_back(nbins);
  param.push_back(ts);
  std::vector<std::vector<double>> p(n);
  std::fill(p.begin(), p.end(), param);
  /*
  std::printf("%-30s %12.8e\n", "tauhat", tauhat);
  std::printf("%-30s %12i\n", "nbins", nbins);
  std::printf("%-30s %12.8e\n", "ts", ts);

    for (std::vector<double>::iterator i = param.begin(); i != param.end(); i++)
    { std::printf("%12.8e\n", *i);

    // print out for debug
    std::for_each(p.begin(), p.end(), [](std::vector<double> &particle) {
      std::printf("%12.8e %12.8e %12.8e \n", particle[0], particle[1],
                  particle[2]);
    });

    std::printf("\n");
  */
  std::transform(data.begin(), data.end(), p.begin(), timecomponent.begin(),
                 ParticleTimesToInteger());
  /*
for (std::vector<int>::iterator i = timecomponent.begin();
  i != timecomponent.end(); i++) {
std::printf("%12i\n", *i);
}
*/
  denseHistogram(timecomponent, histogram, nbins, tauhat);
  /* for (std::vector<int>::iterator i = histogram.begin(); i !=
 histogram.end(); i++) { std::printf("%3i", *i);
   }
 std:
   printf("\n");
 */
  return histogram;
}

/*
 ********************************************************************************
 ********************************************************************************
 * REF:
 * https://stackoverflow.com/questions/14924912/computing-column-sums-of-matrix-vectorvectordouble-with-iterators
 ********************************************************************************
 */
std::vector<double> getColumnMeans(std::vector<std::vector<double>> &dist) {
  std::vector<double> colsums(dist[0].size());

  std::for_each(dist.begin(), dist.end(), [&](const std::vector<double> &row) {
    std::transform(
        row.begin(), row.end(), colsums.begin(), colsums.end(),
        [&dist](double d1, double d2) { return (d1 + d2) / dist.size(); });
  });
  return colsums;
}

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

std::vector<double> vectorSub(std::vector<double> x, std::vector<double> y) {
  std::transform(x.begin(), x.end(), y.begin(), x.begin(),
                 std::minus<double>());
  return x;
}
std::vector<double> CalcRMS(std::vector<std::vector<double>> dist) {
  std::vector<double> avg = getColumnMeans(dist);
  std::vector<std::vector<double>> avgarr;
  std::fill(avgarr.begin(), avgarr.end(), avg);
  std::transform(dist.begin(), dist.end(), avgarr.begin(), dist.begin(),
                 vectorSub);
  std::transform(dist.begin(), dist.end(), dist.begin(), dist.begin(),
                 vectorMultiply);

  std::vector<double> MS = getColumnMeans(dist);
  std::transform(MS.begin(), MS.end(), MS.begin(), (double (*)(double))sqrt);

  return MS;
}
