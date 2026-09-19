// A link/API availability probe, not an EventIO file reader.
#include "initial.h"
#include "io_basic.h"
#include "mc_tel.h"
#include <iostream>

int main() {
  IO_BUFFER* buffer = allocate_io_buffer(4096);
  if (!buffer) return 1;
  // Both public decoding entry points explicitly reject null buffers. Calling
  // them verifies link availability without fabricating an input event.
  int array = 0, telescope = 0, count = 0;
  double photons = 0;
  const int normal = read_tel_photons(nullptr, 0, &array, &telescope, &photons, nullptr, &count);
  const int spatial = read_tel_photons3d(nullptr, 0, &array, &telescope, &photons, nullptr, &count);
  free_io_buffer(buffer);
  if (normal != -1 || spatial != -1) return 2;
  std::cout << "hessio EventIO allocation and 2D/3D photon APIs available\n";
}
