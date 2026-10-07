#pragma once

#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <locale>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace obdeect {

inline bool parse_positive_size(std::string_view text, std::size_t &output) {
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size() && output > 0;
}

inline bool parse_uint32(std::string_view text, std::uint32_t &output) {
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), output);
  return error == std::errc{} && end == text.data() + text.size();
}

// Floating-point from_chars is unavailable in the libc++ shipped with
// Xcode 16. Parse CLI arguments with a fixed locale outside the trace path.
inline bool parse_finite_double(std::string_view text, double &output) {
  std::istringstream input{std::string{text}};
  input.imbue(std::locale::classic());
  input >> std::noskipws >> output;
  return !input.fail() && input.eof() && std::isfinite(output);
}

// A source spectrum is represented by an explicit, equally sampled list of
// wavelengths.  Keeping the parsing at the CLI boundary avoids hidden spectral
// assumptions in the tracing kernel.
inline bool parse_wavelengths_nm(std::string_view text, std::vector<double> &output) {
  if (text.empty() || text.back() == ',')
    return false;
  std::vector<double> parsed;
  std::size_t begin = 0;
  while (begin < text.size()) {
    const std::size_t end = text.find(',', begin);
    const auto token =
        text.substr(begin, end == std::string_view::npos ? text.size() - begin : end - begin);
    double wavelength_nm{};
    if (token.empty() || !parse_finite_double(token, wavelength_nm) || wavelength_nm <= 0.0)
      return false;
    parsed.push_back(wavelength_nm);
    if (end == std::string_view::npos)
      break;
    begin = end + 1;
  }
  if (parsed.empty())
    return false;
  output = std::move(parsed);
  return true;
}

} // namespace obdeect
