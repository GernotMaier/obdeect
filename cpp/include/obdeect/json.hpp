#pragma once

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace obdeect::json {

// Minimal strict JSON document model for optical-model files.
struct Value {
  enum class Kind { null, boolean, number, string, array, object };
  using Array = std::vector<Value>;
  using Object = std::vector<std::pair<std::string, Value>>;

  std::variant<std::monostate, bool, double, std::string, Array, Object> data;

  [[nodiscard]] Kind kind() const { return static_cast<Kind>(data.index()); }
  [[nodiscard]] const Value *find(std::string_view key) const {
    if (kind() != Kind::object)
      return nullptr;
    for (const auto &[name, value] : std::get<Object>(data))
      if (name == key)
        return &value;
    return nullptr;
  }
};

class Parser {
public:
  // Throws std::runtime_error with the byte offset of the first syntax error.
  [[nodiscard]] static Value parse(std::string_view text) {
    Parser parser{text};
    parser.space();
    Value value = parser.value(0);
    parser.space();
    if (parser.at_ != text.size())
      parser.fail("trailing characters");
    return value;
  }

private:
  explicit Parser(std::string_view text) : text_(text) {}

  [[noreturn]] void fail(const char *what) const {
    throw std::runtime_error("JSON syntax error at byte " + std::to_string(at_) + ": " + what);
  }
  void space() {
    while (at_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[at_])))
      ++at_;
  }
  bool take(char c) {
    if (at_ < text_.size() && text_[at_] == c) {
      ++at_;
      return true;
    }
    return false;
  }
  void expect(char c) {
    if (!take(c))
      fail("unexpected character");
  }

  Value value(int depth) {
    if (depth > 64)
      fail("nesting too deep");
    if (at_ == text_.size())
      fail("unexpected end of input");
    switch (text_[at_]) {
    case '{':
      return object(depth);
    case '[':
      return array(depth);
    case '"':
      return Value{string()};
    case 't':
      return literal("true", Value{true});
    case 'f':
      return literal("false", Value{false});
    case 'n':
      return literal("null", Value{});
    default:
      return Value{number()};
    }
  }

  Value literal(std::string_view word, Value result) {
    if (text_.substr(at_, word.size()) != word)
      fail("invalid literal");
    at_ += word.size();
    return result;
  }

  std::string string() {
    expect('"');
    std::string result;
    while (at_ < text_.size()) {
      const char c = text_[at_++];
      if (c == '"')
        return result;
      if (static_cast<unsigned char>(c) < 0x20)
        fail("control character in string");
      if (c != '\\') {
        result.push_back(c);
        continue;
      }
      if (at_ == text_.size())
        break;
      const char e = text_[at_++];
      const std::string_view plain = "\"\\/", escaped = "bfnrt", decoded = "\b\f\n\r\t";
      if (plain.find(e) != std::string_view::npos)
        result.push_back(e);
      else if (const auto i = escaped.find(e); i != std::string_view::npos)
        result.push_back(decoded[i]);
      else
        fail("unsupported escape sequence");
    }
    fail("unterminated string");
  }

  double number() {
    const std::size_t start = at_;
    take('-');
    if (!take('0')) {
      if (at_ == text_.size() || !std::isdigit(static_cast<unsigned char>(text_[at_])))
        fail("invalid number");
      digits();
    }
    if (take('.') && !digits())
      fail("invalid fraction");
    if (take('e') || take('E')) {
      if (!take('+'))
        take('-');
      if (!digits())
        fail("invalid exponent");
    }
    const std::string literal{text_.substr(start, at_ - start)};
    const double result = std::strtod(literal.c_str(), nullptr);
    if (!std::isfinite(result))
      fail("number out of range");
    return result;
  }
  bool digits() {
    const std::size_t start = at_;
    while (at_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[at_])))
      ++at_;
    return at_ > start;
  }

  Value array(int depth) {
    expect('[');
    Value::Array items;
    space();
    if (!take(']')) {
      do {
        space();
        items.push_back(value(depth + 1));
        space();
      } while (take(','));
      expect(']');
    }
    return Value{std::move(items)};
  }

  Value object(int depth) {
    expect('{');
    Value::Object members;
    space();
    if (!take('}')) {
      do {
        space();
        std::string key = string();
        space();
        expect(':');
        space();
        members.emplace_back(std::move(key), value(depth + 1));
        space();
      } while (take(','));
      expect('}');
    }
    return Value{std::move(members)};
  }

  std::string_view text_;
  std::size_t at_{};
};

} // namespace obdeect::json
