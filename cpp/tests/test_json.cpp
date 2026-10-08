#include "obdeect/json.hpp"

#include <cassert>
#include <string>
#include <utility>

int main() {
  using obdeect::json::Value;
  // Bound storage per table entry: inactive JSON containers must not multiply
  // memory consumption for large numeric response tables.
  static_assert(sizeof(Value) <= 64);
  auto root =
      obdeect::json::Parser{
          R"({"null":null,"bool":true,"number":-1.25e+2,"string":"a\n\t\\\"b","array":[0,1.0,2],"object":{}})"}
          .parse();
  assert(root && root->kind() == Value::Kind::object);
  assert(root->find("null")->kind() == Value::Kind::null);
  assert(root->find("bool")->boolean());
  assert(root->find("number")->number() == -125.0);
  assert(root->find("number")->number_text() == "-1.25e+2");
  assert(root->find("string")->string() == "a\n\t\\\"b");
  assert(root->find("array")->array().size() == 3);
  assert(root->find("array")->array()[1].number_text() == "1.0");
  assert(root->find("object")->object().empty());
  assert(!root->find("absent") && !root->find("array")->find("key"));
  auto copy = *root;
  copy.object().clear();
  assert(!root->object().empty());
  auto moved = std::move(*root);
  assert(moved.find("number")->number() == -125.0);
  for (const auto invalid : {"", "[1,]", "{\"a\":}", "01", "1e", "1.", "1e999", "true false",
                             "\"\\u0000\"", "\"unterminated"})
    assert(!obdeect::json::Parser{invalid}.parse());
  assert(obdeect::json::Parser{"[]"}.parse()->array().empty());
  assert(!obdeect::json::Parser{"false"}.parse()->boolean());
  const std::string deep = std::string(130, '[') + "0" + std::string(130, ']');
  assert(!obdeect::json::Parser{deep}.parse());
}
