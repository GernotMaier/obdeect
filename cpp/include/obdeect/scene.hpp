#pragma once

#include "obdeect/artificial_mst.hpp"

namespace obdeect {

enum class TraceMode : std::uint8_t { directed, bounded_nonsequential };

// Immutable compiled-scene seed for the first vertical slice. Future LST/MST
// segmented and SC dual-mirror scene variants share this compile/trace split.
struct CompiledReferenceScene {
  ArtificialMstConfig configuration;
  TraceMode mode{TraceMode::directed};

  [[nodiscard]] bool is_valid() const { return mode == TraceMode::directed && obdeect::is_valid(configuration); }
};

[[nodiscard]] inline std::optional<CompiledReferenceScene> compile_scene(const ArtificialMstConfig& configuration) {
  if (!obdeect::is_valid(configuration)) return std::nullopt;
  return CompiledReferenceScene{configuration};
}

}  // namespace obdeect
