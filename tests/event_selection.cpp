// Test strict score boundaries and exceptional inputs in the production gate.
#define main autofocus_main
#include "../src/cal_ccf.cpp"
#undef main

static void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

int main() try {
  EventSelection selection;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  require(selection.accepts(-1., 1., 0) && selection.accepts(nan, 0., 2),
          "default all mode must preserve the original candidate path");
  selection.selected_only = true;
  for (int component = 0; component < 3; ++component) {
    const double cutoff = selection.minimum_max_mad[component];
    require(!selection.accepts(cutoff, 1., component), "threshold equality must reject");
    require(selection.accepts(std::nextafter(cutoff, inf), 1., component),
            "score immediately above the threshold must pass");
    require(!selection.accepts(std::nextafter(cutoff, 0.), 1., component),
            "score immediately below the threshold must reject");
    require(selection.accepts(cutoff * 2., 1., component) &&
            selection.accepts(cutoff * 2e-20, 1e-20, component),
            "the score must be invariant to a common amplitude scale");
    for (double invalid_mad : {0., -1., nan, inf})
      require(!selection.accepts(100., invalid_mad, component), "invalid MAD must reject");
    for (double invalid_peak : {-1., nan, inf, -inf})
      require(!selection.accepts(invalid_peak, 1., component), "invalid/negative peak must reject");
  }
  selection.minimum_max_mad = {{2., 5., 40.}};
  require(selection.accepts(3., 1., 0) && !selection.accepts(3., 1., 1),
          "R and T thresholds must be independent");
  require(!selection.accepts(39., 1., 2) && selection.accepts(41., 1., 2),
          "custom U threshold");
  std::cout << "PASS: selection defaults, strict boundaries, scale, finite inputs and independent thresholds\n";
} catch (const std::exception &error) {
  std::cerr << error.what() << '\n';
  return 1;
}
