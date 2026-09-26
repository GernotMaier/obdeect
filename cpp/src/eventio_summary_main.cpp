#include "obdeect/eventio_photon_input.hpp"

#include <array>
#include <cstdint>
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: obdeect_eventio_summary CORSIKA_TELFIL\n";
    return 2;
  }
  try {
    obdeect::EventioPhotonReader reader(argv[1]);
    std::array<obdeect::OpticalPhoton, 65536> bunches{};
    std::uint64_t batches = 0, count = 0;
    double photons = 0;
    for (;;) {
      const auto result = reader.read(bunches);
      if (result.count > 0) {
        ++batches;
        count += result.count;
        for (std::size_t index = 0; index < result.count; ++index)
          photons += bunches[index].weight;
      }
      if (result.eof) break;
    }
    std::cout << "CORSIKA EventIO batches=" << batches << " bunches=" << count
              << " photons=" << photons << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
