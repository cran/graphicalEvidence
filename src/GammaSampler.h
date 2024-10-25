#pragma once

/* Class for generating a gamma sample with input shape and scale */


class GammaSampler {
public:
  GammaSampler()
  {
    std::mt19937 gen(rd());
  }

  void SetSeed(unsigned int seed) {
    gen.seed(seed);
  }

  double GetSample(const double shape, const double scale) {
    std::gamma_distribution<double> distribution(shape, scale);
    return distribution(gen);
  }

private:
  std::random_device rd;
  std::mt19937 gen;
};

/* Global gamma sampler */
extern GammaSampler g_rgamma;
