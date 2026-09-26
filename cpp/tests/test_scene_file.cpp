#include "obdeect/scene_file.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>

int main() {
  const char* path = "obdeect-native-scene-test.csv";
  {
    std::ofstream output(path);
    output << "obdeect-scene-v1\n"
              "provenance,LSTN-design,7.0.0,aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n"
              "surface_id,role,shape,cx_m,cy_m,cz_m,nx,ny,nz,tx,ty,tz,diameter_m\n"
              "0,mirror,circle,0,0,0,0,0,1,1,0,0,1\n"
              "1,detector,circle,0,0,10,0,0,1,1,0,0,1\n";
  }
  const auto scene = obdeect::read_native_scene(path);
  assert(scene.has_value());
  assert(scene->primary_facets.size() == 1);
  assert(scene->detector_surfaces.size() == 1);
  std::remove(path);
}
