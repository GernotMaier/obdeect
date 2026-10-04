#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace obdeect::detail {

// Isolate real roots on a finite interval using derivative roots. Between two
// consecutive stationary points a polynomial is monotone, so bisection cannot
// skip a crossing. Stationary roots also retain even-multiplicity tangencies.
// Fixed storage bounds the degree and keeps the transport path allocation-free.
using Polynomial = std::array<long double, 27>;
struct PolynomialRoots {
  std::array<long double, 26> values{};
  std::size_t count{};
};

[[nodiscard]] inline long double polynomial_value(const Polynomial &coefficients,
                                                  std::size_t degree, long double x) {
  long double value = coefficients[degree];
  while (degree != 0)
    value = value * x + coefficients[--degree];
  return value;
}

[[nodiscard]] inline PolynomialRoots polynomial_roots_unit_interval(Polynomial coefficients,
                                                                    std::size_t degree) {
  while (degree != 0 && coefficients[degree] == 0)
    --degree;
  PolynomialRoots result{};
  if (degree == 0)
    return result;
  long double scale = 0;
  for (std::size_t i = 0; i <= degree; ++i)
    scale = std::max(scale, std::abs(coefficients[i]));
  if (!std::isfinite(scale) || scale == 0)
    return result;
  for (std::size_t i = 0; i <= degree; ++i)
    coefficients[i] /= scale;
  if (degree == 1) {
    const long double root = -coefficients[0] / coefficients[1];
    if (root >= 0 && root <= 1)
      result.values[result.count++] = root;
    return result;
  }
  Polynomial derivative{};
  for (std::size_t i = 1; i <= degree; ++i)
    derivative[i - 1] = coefficients[i] * i;
  const auto stationary = polynomial_roots_unit_interval(derivative, degree - 1);
  constexpr long double tolerance = 128 * std::numeric_limits<long double>::epsilon();
  const auto add = [&](long double root) {
    if (result.count < result.values.size() &&
        (result.count == 0 || std::abs(root - result.values[result.count - 1]) > tolerance))
      result.values[result.count++] = root;
  };
  long double left = 0;
  long double left_value = polynomial_value(coefficients, degree, left);
  if (std::abs(left_value) <= tolerance)
    add(left);
  for (std::size_t i = 0; i <= stationary.count; ++i) {
    const long double right = i == stationary.count ? 1 : stationary.values[i];
    const long double right_value = polynomial_value(coefficients, degree, right);
    if (std::abs(left_value) > tolerance && std::abs(right_value) > tolerance &&
        std::signbit(left_value) != std::signbit(right_value)) {
      long double low = left, high = right, low_value = left_value;
      for (int iteration = 0; iteration < 100; ++iteration) {
        const long double middle = (low + high) / 2;
        if (middle == low || middle == high)
          break;
        const long double value = polynomial_value(coefficients, degree, middle);
        if (value == 0) {
          low = high = middle;
          break;
        }
        if (std::signbit(value) == std::signbit(low_value)) {
          low = middle;
          low_value = value;
        } else
          high = middle;
      }
      add((low + high) / 2);
    }
    if (std::abs(right_value) <= tolerance)
      add(right);
    left = right;
    left_value = right_value;
  }
  return result;
}

} // namespace obdeect::detail
