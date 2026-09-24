#include "obdeect/photon_input.hpp"

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(1);
  }
}

class TemporaryFile {
 public:
  explicit TemporaryFile(const std::string& contents)
      : path_(std::filesystem::temp_directory_path() / "obdeect-photon-input-test.csv") {
    std::ofstream output(path_);
    output << contents;
  }
  ~TemporaryFile() { std::filesystem::remove(path_); }
  [[nodiscard]] std::string path() const { return path_.string(); }

 private:
  std::filesystem::path path_;
};

constexpr const char* kHeader =
    "run_id,event_id,array_id,telescope_id,photon_id,x_m,y_m,z_m,dx,dy,dz,wavelength_nm,time_ns,weight,"
    "bunch_id,emission_height_m,emission_distance_m,telescope_x_m,telescope_y_m,telescope_z_m,array_reuse_weight\n";

}  // namespace

int main() {
  using namespace obdeect;

  // T-INP-001: all raw-bunch and event/telescope metadata survives CSV
  // ingestion, including lambda=0 as an unresolved spectrum sentinel.
  TemporaryFile input{std::string{kHeader} +
                      "7,8,9,10,101,1,2,3,0,0,-1,0,4.5,0.25,500,12000,3000,4,5,6,1.5\n" +
                      "7,8,9,10,102,2,3,4,0,0,-1,350,5.5,0.75,501,13000,3100,4,5,6,1.5\n" +
                      "7,11,9,10,103,3,4,5,0,0,-1,400,6.5,1.0,502,14000,3200,4,5,6,2.0\n"};
  CsvPhotonReader reader(input.path());
  std::array<OpticalPhoton, 1> batch{};
  const auto first = reader.read(batch);
  require(first.count == 1 && !first.eof, "capacity chunk is not falsely EOF");
  require(first.context.run_id == 7 && first.context.event_id == 8 && first.context.telescope_id == 10,
          "event/telescope identifiers preserved");
  require(first.context.telescope_position_m.x == 4.0 && first.context.telescope_position_m.y == 5.0 &&
              first.context.telescope_position_m.z == 6.0 && first.context.array_reuse_weight == 1.5,
          "telescope position and fractional reuse weight preserved");
  require(batch[0].photon_id == 101 && batch[0].wavelength_nm == 0.0 && batch[0].bunch_id == 500 &&
              batch[0].emission_height_m == 12000.0 && batch[0].emission_distance_m == 3000.0,
          "raw bunch provenance and unspecified wavelength preserved");

  const auto second = reader.read(batch);
  require(second.count == 1 && !second.eof && batch[0].photon_id == 102,
          "pending row preserves order across capacity chunks");
  const auto third = reader.read(batch);
  require(third.count == 1 && third.eof && third.context.event_id == 11 && batch[0].photon_id == 103,
          "event boundary begins a new batch and final row is retained");
  require(reader.read(batch).count == 0 && reader.read(batch).eof, "repeated EOF is stable");

  // T-INP-003: CRLF records are accepted without retaining the carriage
  // return in the final header or numeric field.
  std::string crlf_header{kHeader};
  for (std::size_t position = 0; (position = crlf_header.find('\n', position)) != std::string::npos;
       position += 2) {
    crlf_header.insert(position, "\r");
  }
  TemporaryFile crlf_input{crlf_header +
                           "7,8,9,10,101,1,2,3,0,0,-1,0,4.5,0.25,500,12000,3000,4,5,6,1.5\r\n"};
  CsvPhotonReader crlf_reader(crlf_input.path());
  const auto crlf_result = crlf_reader.read(batch);
  require(crlf_result.count == 1 && crlf_result.eof && batch[0].photon_id == 101 &&
              batch[0].weight == 0.25,
          "CRLF CSV records are parsed without carriage-return suffixes");

  // T-INP-002: required schema and invalid values fail explicitly rather than
  // becoming empty events or silently repaired photons.
  TemporaryFile missing{"run_id,event_id\n1,2\n"};
  bool rejected = false;
  try {
    CsvPhotonReader invalid(missing.path());
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "missing required column is rejected");

  TemporaryFile malformed{std::string{kHeader} +
                          "1,2,3,4,5,0,0,0,0,0,0,400,0,1,0,0,0,0,0,0,1\n"};
  rejected = false;
  try {
    CsvPhotonReader invalid(malformed.path());
    (void)invalid.read(batch);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "zero direction is rejected");

  TemporaryFile negative_id{std::string{kHeader} +
                            "1,2,3,4,-5,0,0,0,0,0,-1,400,0,1,0,0,0,0,0,0,1\n"};
  rejected = false;
  try {
    CsvPhotonReader invalid(negative_id.path());
    (void)invalid.read(batch);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "negative unsigned identifiers are rejected rather than wrapped");

  std::cout << "photon input tests passed\n";
}
