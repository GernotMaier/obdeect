#pragma once

#include <charconv>
#include <cmath>
#include <cstddef>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>

namespace obdeect {

inline bool parse_positive_size(std::string_view text, std::size_t& output) {
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size() && output > 0;
}

// Floating-point from_chars is unavailable in the libc++ shipped with
// Xcode 16. Parse CLI arguments with a fixed locale outside the trace path.
inline bool parse_finite_double(std::string_view text, double& output) {
  std::istringstream input{std::string{text}};
  input.imbue(std::locale::classic());
  input >> std::noskipws >> output;
  return !input.fail() && input.eof() && std::isfinite(output);
}

}  // namespace obdeect
