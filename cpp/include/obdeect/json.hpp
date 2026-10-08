#pragma once

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::json {

// Small, strict JSON reader for versioned optical-model files.  Keeping this
// boundary standard-library-only avoids imposing a JSON dependency on the
// tracing kernels or on binary-wheel users.
struct Value {
  enum class Kind { null, boolean, number, string, array, object };
  struct Number {
    double value;
    // Preserve number spelling for canonical artifact hashing.
    std::string text;
  };
  using Array = std::vector<Value>;
  using Object = std::vector<std::pair<std::string, Value>>;

  Value() = default;
  explicit Value(bool value) : data_(value) {}
  explicit Value(Number value) : data_(std::move(value)) {}
  explicit Value(std::string value) : data_(std::move(value)) {}
  explicit Value(Array value) : data_(std::move(value)) {}
  explicit Value(Object value) : data_(std::move(value)) {}

  [[nodiscard]] Kind kind() const { return static_cast<Kind>(data_.index()); }
  [[nodiscard]] bool boolean() const { return std::get<bool>(data_); }
  bool &boolean() { return std::get<bool>(data_); }
  [[nodiscard]] double number() const { return std::get<Number>(data_).value; }
  [[nodiscard]] const std::string &number_text() const { return std::get<Number>(data_).text; }
  [[nodiscard]] const std::string &string() const { return std::get<std::string>(data_); }
  std::string &string() { return std::get<std::string>(data_); }
  [[nodiscard]] const Array &array() const { return std::get<Array>(data_); }
  Array &array() { return std::get<Array>(data_); }
  [[nodiscard]] const Object &object() const { return std::get<Object>(data_); }
  Object &object() { return std::get<Object>(data_); }

  [[nodiscard]] const Value *find(std::string_view key) const {
    if (kind() != Kind::object)
      return nullptr;
    for (const auto &[name, value] : object())
      if (name == key)
        return &value;
    return nullptr;
  }

private:
  // Only the active JSON type occupies storage; numeric table entries do not
  // carry unused string, array and object containers.
  std::variant<std::monostate, bool, Number, std::string, Array, Object> data_;
};

class Parser {
public:
  explicit Parser(std::string_view input) : input_(input) {}

  [[nodiscard]] std::optional<Value> parse() {
    skip_space();
    auto value = parse_value(0);
    skip_space();
    return value && position_ == input_.size() ? value : std::nullopt;
  }

private:
  static constexpr std::size_t kMaximumDepth = 128;

  void skip_space() {
    while (position_ < input_.size() && std::isspace(static_cast<unsigned char>(input_[position_])))
      ++position_;
  }

  [[nodiscard]] bool consume(char character) {
    if (position_ == input_.size() || input_[position_] != character)
      return false;
    ++position_;
    return true;
  }

  [[nodiscard]] std::optional<Value> parse_value(std::size_t depth) {
    if (depth > kMaximumDepth || position_ == input_.size())
      return std::nullopt;
    switch (input_[position_]) {
    case '{':
      return parse_object(depth + 1);
    case '[':
      return parse_array(depth + 1);
    case '"': {
      auto text = parse_string();
      if (!text)
        return std::nullopt;
      return Value{std::move(*text)};
    }
    case 't':
      return consume_literal("true", Value::Kind::boolean, true);
    case 'f':
      return consume_literal("false", Value::Kind::boolean, false);
    case 'n':
      return consume_literal("null", Value::Kind::null);
    default:
      return parse_number();
    }
  }

  [[nodiscard]] std::optional<Value> consume_literal(std::string_view literal, Value::Kind kind,
                                                     bool boolean = false) {
    if (input_.substr(position_, literal.size()) != literal)
      return std::nullopt;
    position_ += literal.size();
    if (kind == Value::Kind::boolean)
      return Value{boolean};
    return Value{};
  }

  [[nodiscard]] std::optional<std::string> parse_string() {
    if (!consume('"'))
      return std::nullopt;
    std::string result;
    while (position_ < input_.size()) {
      const char character = input_[position_++];
      if (character == '"')
        return result;
      if (static_cast<unsigned char>(character) < 0x20)
        return std::nullopt;
      if (character != '\\') {
        result.push_back(character);
        continue;
      }
      if (position_ == input_.size())
        return std::nullopt;
      switch (input_[position_++]) {
      case '"':
        result.push_back('"');
        break;
      case '\\':
        result.push_back('\\');
        break;
      case '/':
        result.push_back('/');
        break;
      case 'b':
        result.push_back('\b');
        break;
      case 'f':
        result.push_back('\f');
        break;
      case 'n':
        result.push_back('\n');
        break;
      case 'r':
        result.push_back('\r');
        break;
      case 't':
        result.push_back('\t');
        break;
      // Model schema keys and values are UTF-8.  Escaped Unicode is not
      // emitted by the compiler and is rejected rather than mis-decoded.
      default:
        return std::nullopt;
      }
    }
    return std::nullopt;
  }

  [[nodiscard]] std::optional<Value> parse_number() {
    const std::size_t start = position_;
    if (position_ < input_.size() && input_[position_] == '-')
      ++position_;
    if (position_ == input_.size())
      return std::nullopt;
    if (input_[position_] == '0') {
      ++position_;
    } else if (input_[position_] >= '1' && input_[position_] <= '9') {
      do {
        ++position_;
      } while (position_ < input_.size() &&
               std::isdigit(static_cast<unsigned char>(input_[position_])));
    } else {
      return std::nullopt;
    }
    if (position_ < input_.size() && input_[position_] == '.') {
      ++position_;
      const std::size_t fraction = position_;
      while (position_ < input_.size() &&
             std::isdigit(static_cast<unsigned char>(input_[position_])))
        ++position_;
      if (fraction == position_)
        return std::nullopt;
    }
    if (position_ < input_.size() && (input_[position_] == 'e' || input_[position_] == 'E')) {
      ++position_;
      if (position_ < input_.size() && (input_[position_] == '+' || input_[position_] == '-'))
        ++position_;
      const std::size_t exponent = position_;
      while (position_ < input_.size() &&
             std::isdigit(static_cast<unsigned char>(input_[position_])))
        ++position_;
      if (exponent == position_)
        return std::nullopt;
    }
    const std::string text{input_.substr(start, position_ - start)};
    char *end = nullptr;
    const double value = std::strtod(text.c_str(), &end);
    if (end != text.c_str() + text.size() || !std::isfinite(value))
      return std::nullopt;
    return Value{Value::Number{value, text}};
  }

  [[nodiscard]] std::optional<Value> parse_array(std::size_t depth) {
    if (!consume('['))
      return std::nullopt;
    Value result{Value::Array{}};
    skip_space();
    if (consume(']'))
      return result;
    while (true) {
      skip_space();
      auto value = parse_value(depth);
      if (!value)
        return std::nullopt;
      result.array().push_back(std::move(*value));
      skip_space();
      if (consume(']'))
        return result;
      if (!consume(','))
        return std::nullopt;
    }
  }

  [[nodiscard]] std::optional<Value> parse_object(std::size_t depth) {
    if (!consume('{'))
      return std::nullopt;
    Value result{Value::Object{}};
    skip_space();
    if (consume('}'))
      return result;
    while (true) {
      skip_space();
      auto key = parse_string();
      if (!key)
        return std::nullopt;
      skip_space();
      if (!consume(':'))
        return std::nullopt;
      skip_space();
      auto value = parse_value(depth);
      if (!value)
        return std::nullopt;
      result.object().emplace_back(std::move(*key), std::move(*value));
      skip_space();
      if (consume('}'))
        return result;
      if (!consume(','))
        return std::nullopt;
    }
  }

  std::string_view input_;
  std::size_t position_{};
};

} // namespace obdeect::json
