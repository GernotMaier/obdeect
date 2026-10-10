"""Reviewed non-transport parameter roles for simulation-models records.

The optical endpoint is the physical sensor surface. Conversion to electrons,
source generation, digitization and triggers retain their upstream/downstream
ownership. Names absent from this explicit catalogue remain unsupported.
"""

NON_TRANSPORT_PARAMETER_ROLES = {}

for _name in (
    "calibration_devices",
    "nsb_autoscale_airmass",
    "nsb_gain_drop_scale",
    "nsb_offaxis",
    "nsb_pixel_rate",
):
    NON_TRANSPORT_PARAMETER_ROLES[_name] = "source generation or calibration-device inventory"

for _name in (
    "pixel_cells",
    "pm_collection_efficiency",
    "pm_photoelectron_spectrum",
    "qe_variation",
    "quantum_efficiency",
):
    NON_TRANSPORT_PARAMETER_ROLES[_name] = (
        "photon-to-electron conversion after the optical endpoint"
    )

for _name in ("array_element_position_utm",):
    NON_TRANSPORT_PARAMETER_ROLES[_name] = (
        "alternate geodetic coordinates; ground coordinates define the optical input frame"
    )

for _name in (
    "adjust_gain",
    "asum_clipping",
    "asum_offset",
    "asum_shaping",
    "asum_threshold",
    "camera_trigger_groups",
    "camera_trigger_members",
    "channels_per_chip",
    "default_trigger",
    "disc_ac_coupled",
    "disc_bins",
    "disc_start",
    "discriminator_amplitude",
    "discriminator_fall_time",
    "discriminator_gate_length",
    "discriminator_hysteresis",
    "discriminator_output_amplitude",
    "discriminator_output_var_percent",
    "discriminator_pulse_shape",
    "discriminator_rise_time",
    "discriminator_scale_threshold",
    "discriminator_sigsum_over_threshold",
    "discriminator_threshold",
    "discriminator_time_over_threshold",
    "discriminator_var_gate_length",
    "discriminator_var_sigsum_over_threshold",
    "discriminator_var_threshold",
    "discriminator_var_time_over_threshold",
    "dsum_clipping",
    "dsum_ignore_below",
    "dsum_offset",
    "dsum_pedsub",
    "dsum_pre_clipping",
    "dsum_prescale",
    "dsum_presum_max",
    "dsum_presum_shift",
    "dsum_shaping",
    "dsum_shaping_renormalize",
    "dsum_threshold",
    "dsum_zero_clip",
    "fadc_ac_coupled",
    "fadc_amplitude",
    "fadc_bins",
    "fadc_compensate_pedestal",
    "fadc_dev_pedestal",
    "fadc_err_compensate_pedestal",
    "fadc_err_pedestal",
    "fadc_lg_amplitude",
    "fadc_lg_compensate_pedestal",
    "fadc_lg_dev_pedestal",
    "fadc_lg_err_compensate_pedestal",
    "fadc_lg_err_pedestal",
    "fadc_lg_max_signal",
    "fadc_lg_max_sum",
    "fadc_lg_noise",
    "fadc_lg_pedestal",
    "fadc_lg_sensitivity",
    "fadc_lg_sysvar_pedestal",
    "fadc_lg_var_pedestal",
    "fadc_lg_var_sensitivity",
    "fadc_long_event_threshold",
    "fadc_long_sum_bins",
    "fadc_long_sum_offset",
    "fadc_max_signal",
    "fadc_max_sum",
    "fadc_mhz",
    "fadc_noise",
    "fadc_pedestal",
    "fadc_pulse_shape",
    "fadc_sensitivity",
    "fadc_sum_bins",
    "fadc_sum_offset",
    "fadc_sysvar_pedestal",
    "fadc_var_pedestal",
    "fadc_var_sensitivity",
    "flatfielding",
    "gain_variation",
    "hg_lg_variation",
    "multiplicity_offset",
    "muon_mono_threshold",
    "num_gains",
    "only_triggered_telescopes",
    "photon_delay",
    "pixeltrg_time_step",
    "pm_average_gain",
    "pm_gain_index",
    "pm_transit_time",
    "pm_voltage_variation",
    "random_mono_probability",
    "teltrig_min_sigsum",
    "teltrig_min_time",
    "transit_time_calib_error",
    "transit_time_compensate_error",
    "transit_time_compensate_step",
    "transit_time_error",
    "transit_time_jitter",
    "transit_time_random",
    "trigger_current_limit",
    "trigger_delay_compensation",
    "trigger_pixels",
):
    NON_TRANSPORT_PARAMETER_ROLES[_name] = (
        "sensor/electronics calibration, digitization or trigger after optical transport"
    )

del _name
