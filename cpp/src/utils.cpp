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
    // std::printf("%-30s %12.8e\n", "dtsamp2", dtsamp2);
    out = (data[4] - params[2] + params[0]) / dtsamp2;
    // std::printf("%-30s %12.8e\n", "data[4]", data[4]);
    // std::printf("%-30s %12.8e\n", "out doub", out);
    out = (int)(out + 0.5);
    // std::printf("%-30s %12.8e\n", "out int", out);
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

  for (std::vector<double>::iterator i = param.begin(); i != param.end(); i++) {
    std::printf("%12.8e\n", *i);
  }

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
  /*
  for (std::vector<int>::iterator i = histogram.begin(); i != histogram.end();
       i++) {
    std::printf("%3i", *i);
  };
  std::printf("\n");
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
        row.begin(), row.end(), colsums.begin(), colsums.begin(),
        [&](double d1, double d2) { return (d1 + d2) / dist.size(); });
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
std::vector<double> CalcRMS(std::vector<std::vector<double>> &dist) {
  std::vector<double> avg = getColumnMeans(dist);

  // for (std::vector<double>::const_iterator i = avg.begin(); i != avg.end();
  // ++i)
  //  std::printf("%-30s %12.8e\n", "avg", *i);

  std::vector<std::vector<double>> avgarr(dist.size());
  std::fill(avgarr.begin(), avgarr.end(), avg);

  // std::for_each(avgarr.begin(), avgarr.end(), [](std::vector<double> &a) {
  //  std::printf("%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", a[0], a[1],
  //  a[2],
  //              a[3], a[4], a[5]);
  //});
  std::vector<std::vector<double>> distcopy = dist;
  std::transform(dist.begin(), dist.end(), avgarr.begin(), distcopy.begin(),
                 vectorSub);
  std::transform(distcopy.begin(), distcopy.end(), distcopy.begin(),
                 distcopy.begin(), vectorMultiply);

  std::vector<double> MS = getColumnMeans(distcopy);
  std::transform(MS.begin(), MS.end(), MS.begin(), (double (*)(double))sqrt);

  return MS;
}

std::vector<double> HistogramToSQRTofCumul(std::vector<int> inputHistogram,
                                           double coeff) {

  /* The name is a badly chose here as it was originally used to calculate
   * cumulated distributions which      */
  /* turned out to be not necessary - the function takes a histogram as input,
   * multiplies it with a constant  */
  /* vector before taking the sqrt of each element. This produces a vector
   * that is used in the IBSNew routine */
  /* to multiply with particle momenta representing the IBS contribution */
  int n = inputHistogram.size();

  std::vector<double> vcoeff(n);
  std::fill(vcoeff.begin(), vcoeff.end(), coeff);
  /*
  for (std::vector<double>::const_iterator i = vcoeff.begin();
       i != vcoeff.end(); ++i)
    std::printf("%-30s %12.8e\n", "avg", *i);
    */
  // fill constant vector

  // multiply with constant
  std::transform(inputHistogram.begin(), inputHistogram.end(), vcoeff.begin(),
                 vcoeff.begin(), std::multiplies<double>());

  // take sqrt
  std::transform(vcoeff.begin(), vcoeff.end(), vcoeff.begin(),
                 (double (*)(double))sqrt);

  return vcoeff;
}

std::map<std::string, double> readInput(std::string filename) {
  std::vector<std::string> ALLOWEDKEYS{
      "bucket",      "atomNumber", "charge",    "nMacro",
      "nReal",       "nbins",      "sigs",      "seed",
      "ex",          "ey",         "timeRatio", "fracibstot",
      "ibsCoupling", "model",      "nturns",    "nwrite"};

  std::map<std::string, double> out;
  std::string line;
  std::ifstream file(filename);

  std::getline(file, line);

  // check if file is open
  if (file.is_open()) {
    std::string key;
    double value;
    std::istringstream iss(line);
    iss >> key >> value;

    vector<string>::iterator it =
        find(ALLOWEDKEYS.begin(), ALLOWEDKEYS.end(), key);
    // cout << key << " " << value << " " << endl;
    if (it != ALLOWEDKEYS.end()) {
      out[key] = value;
    }
    while (!file.eof()) {
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> key >> value;

      vector<string>::iterator it =
          find(ALLOWEDKEYS.begin(), ALLOWEDKEYS.end(), key);
      // cout << key << " " << value << " " << endl;
      if (it != ALLOWEDKEYS.end()) {
        out[key] = value;
      }
    }

    file.close();
  }
  return out;
}

bool isInLong(std::vector<double> particle, double tauhat, double synctime) {
  return !(abs(particle[4] - synctime) < tauhat);
}