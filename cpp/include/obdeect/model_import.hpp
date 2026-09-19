#pragma once

#include "obdeect/ctao_models.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace obdeect {

// Canonical import target. Parsers live outside the hot C++ kernels and must
// populate provenance before a model can be compiled into a scene.
struct ModelProvenance {
  std::string model_name;
  std::string model_version;
  std::string content_hash;
};

struct ImportedOpticalModel {
  TelescopeOpticalFamily family;
  ModelProvenance provenance;
  bool requires_segment_list{true};
  bool requires_camera_model{true};
};

// Import boundary for a resolved CTAO model selection. A file-format adapter
// (JSON, a database export, etc.) owns parsing and supplies immutable
// provenance; tracing never discovers an arbitrary model implicitly.
struct ImportedCtaoReferenceModel {
  CtaoReferenceModel optical;
  ImportedOpticalModel import;
};

[[nodiscard]] inline std::optional<ImportedCtaoReferenceModel> import_ctao_reference_model(
    std::string_view name, ModelProvenance provenance) {
  const auto optical = ctao_reference_model(name);
  if (!optical || provenance.model_name.empty() || provenance.model_version.empty() ||
      provenance.content_hash.size() != 64 || provenance.model_name != name) {
    return std::nullopt;
  }
  for (const char character : provenance.content_hash) {
    if (!((character >= '0' && character <= '9') || (character >= 'a' && character <= 'f'))) {
      return std::nullopt;
    }
  }
  return ImportedCtaoReferenceModel{*optical, {optical->family, provenance, true, true}};
}

}  // namespace obdeect
