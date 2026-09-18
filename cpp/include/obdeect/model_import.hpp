#pragma once

#include "obdeect/ctao_optical_specs.hpp"

#include <string_view>

namespace obdeect {

// Canonical import target. Parsers live outside the hot C++ kernels and must
// populate provenance before a model can be compiled into a scene.
struct ModelProvenance {
  std::string_view model_name;
  std::string_view model_version;
  std::string_view content_hash;
};

struct ImportedOpticalModel {
  TelescopeOpticalFamily family;
  ModelProvenance provenance;
  bool requires_segment_list{true};
  bool requires_camera_model{true};
};

}  // namespace obdeect
