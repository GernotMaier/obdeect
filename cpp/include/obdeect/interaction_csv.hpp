#pragma once

#include "obdeect/photon_buffer.hpp"
#include "obdeect/photon_input.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace obdeect {

// Optional diagnostics are boundary I/O. Both caps include all stored records;
// the byte cap includes the header. Transport never allocates or writes here.
class InteractionCsvWriter {
public:
  InteractionCsvWriter(const std::string &path, std::size_t record_limit, std::size_t byte_limit)
      : output_(path), record_limit_(record_limit), byte_limit_(byte_limit) {
    const std::string header =
        "contract_version,frame,run_id,event_id,array_id,telescope_id,photon_id,bunch_id,"
        "interaction_index,surface_id,kind,x_m,y_m,z_m,nx,ny,nz,in_dx,in_dy,in_dz,"
        "out_dx,out_dy,out_dz,segment_length_m,geometric_path_m,optical_path_m,"
        "time_ns,wavelength_nm,incident_weight,outgoing_weight,terminal_status\n";
    if (!output_ || record_limit == 0 || byte_limit < header.size())
      throw std::invalid_argument(
          "cannot open interaction diagnostics or byte cap is smaller than header");
    output_ << header;
    bytes_written_ = header.size();
  }

  template <std::size_t Capacity>
  void write(const PhotonBatchContext &context, const OpticalPhoton &photon,
             const BasicPathRecord<Capacity> &path) {
    double cumulative = 0, previous_throughput = 1;
    for (std::size_t i = 0; i + 1 < path.point_count; ++i) {
      const double segment = norm(path.points_m[i + 1] - path.points_m[i]);
      cumulative += segment;
      if (records_written_ == record_limit_) {
        truncated_ = true;
        return;
      }
      const auto &point = path.points_m[i + 1];
      const auto &normal = path.interaction_normals[i];
      const auto &incoming = path.interaction_incoming_directions[i];
      const auto &outgoing = path.interaction_outgoing_directions[i];
      const auto kind = path.interaction_kinds[i];
      const char *name = kind == OpticalInteractionKind::mirror     ? "mirror"
                         : kind == OpticalInteractionKind::detector ? "detector"
                         : kind == OpticalInteractionKind::obscurer ? "obscurer"
                                                                    : "refractive_interface";
      std::ostringstream record;
      record << std::setprecision(17) << "obdeect-interaction-v1,telescope_local," << context.run_id
             << ',' << context.event_id << ',' << context.array_id << ',' << context.telescope_id
             << ',' << photon.photon_id << ',' << photon.bunch_id << ',' << i << ','
             << path.interaction_surface_ids[i] << ',' << name << ',' << point.x << ',' << point.y
             << ',' << point.z << ',' << normal.x << ',' << normal.y << ',' << normal.z << ','
             << incoming.x << ',' << incoming.y << ',' << incoming.z << ',' << outgoing.x << ','
             << outgoing.y << ',' << outgoing.z << ',' << segment << ',' << cumulative << ','
             << (path.material_transport ? path.interaction_optical_path_m[i] : cumulative) << ','
             << photon.time_ns + (path.material_transport ? path.interaction_group_delay_ns[i]
                                                          : cumulative / kSpeedOfLightMPerNs)
             << ',' << photon.wavelength_nm << ',' << photon.weight * previous_throughput << ','
             << (kind == OpticalInteractionKind::obscurer
                     ? 0
                     : photon.weight * path.interaction_throughput[i])
             << ',' << to_string(path.status) << '\n';
      const auto text = record.str();
      if (text.size() > byte_limit_ - bytes_written_) {
        truncated_ = true;
        return;
      }
      output_ << text;
      if (!output_)
        throw std::runtime_error("failed writing interaction diagnostics");
      bytes_written_ += text.size();
      ++records_written_;
      previous_throughput = path.interaction_throughput[i];
    }
  }

  [[nodiscard]] std::size_t records_written() const { return records_written_; }
  [[nodiscard]] bool truncated() const { return truncated_; }

private:
  std::ofstream output_;
  std::size_t record_limit_, byte_limit_, records_written_{}, bytes_written_{};
  bool truncated_{};
};

} // namespace obdeect
