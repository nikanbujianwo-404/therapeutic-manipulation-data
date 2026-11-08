# Data dictionary

The executable workflow uses the frozen inputs described here. Numeric values were not changed to match manuscript results. `schemas/FILE_SCHEMAS.json` lists the actual fields, dimensions/dtypes and row counts of every public data file. `schemas/FIELD_DICTIONARY.csv` provides a field-level catalog; the definitions below describe the scientific relationships between fields. An unknown vendor convention is explicitly marked unknown rather than inferred from a device's usual settings.

The `source` column in result tables retains analysis-relative paths such as `evidence/formal_ati_imu_phase_profiles.csv`: resolve these relative to `data/processed/`, or to the new reproduction output directory, not to the repository root. Manifest `*_path` fields are instead repository-root relative.

## Records and layout

`ExpertA_P3_B1_R2` means practitioner ExpertA, specimen P3, round 1, repetition 2. `record_manifest.csv` contains one row for every cell in the 2 × 4 × 2 × 2 design. `practitioner`, `specimen`, `round`, and `repetition` encode those factors. Original analysis tables retain the column names `expert`, `material`, `block`, and `repeat`, with the same meanings and public values. `record_a`/`record_b` identify an ordered pair; `pair_id` joins them. `*_path` columns are relative to the repository root after extraction; shared tables are selected using `record_id` or the two pair members. `status` and `notes` describe availability, not measurement validity.

P1 is Ecoflex 00-10; P3 is Ecoflex 00-20; P2 is Ecoflex 00-30; P4 is Dragon Skin FX-Pro. The measured mechanical order is P1 < P3 < P2 < P4. `B1/B2` are rounds, not independently verified remounting sessions. `R1/R2` are repetitions within each practitioner–specimen–round cell. `E1/E2` in the example figures designate selected executions, whose record IDs are retained in figure-source tables.

`comparison=expert`: 16 matched practitioner pairs. `repeat`: 16 same-condition repetition pairs. `material`: 48 specimen pairs within the same practitioner, round and repetition. A pair contributes nine scalar channels; 720 scalar-pair rows represent 80 pairs, not 720 independent participants. The 32 recordings are nested within two practitioners.

Material trial IDs, such as `P1_indentation_R1`, identify a separate 4 × 3 experiment; `material_manifest.csv` indexes its 12 trials, source paths, and the original `max_sync_gap_sec=0.02` used for measured-force interpolation.

## Raw acquisition files

### ATI Mini45

`raw/records/{record_id}/ati/decoded.csv` is a byte-identical copy of the full decoded acquisition table. Output rate is **200 Hz**, as recorded in `sensor_fields.comrdtrate`. The **7000 Hz** `runrate` is the internal force/torque counter basis, not the output sampling rate.

| Columns | Meaning / units |
|---|---|
| `host_receive_monotonic_ns`, `host_receive_unix_ns` | Software receipt timestamps, respectively monotonic and Unix nanoseconds |
| `host_elapsed_sec` | Receipt time in seconds relative to the recording's common software start gate |
| `rdt_sequence`, `ft_sequence` | Packet and internal force/torque sequence counters; do not interchange their tick rates |
| `status`, `sequence_delta`, `valid`, `bias_id` | Device status word, packet sequence increment, acquisition validity flag (0/1), applied tare index |
| `raw_fx_count`…`raw_fz_count`, `raw_tx_count`…`raw_tz_count` | Signed digital counts in Fx, Fy, Fz, Tx, Ty, Tz order |
| `fx_N`, `fy_N`, `fz_N` | Count-scaled forces in N |
| `tx_Nm`, `ty_Nm`, `tz_Nm` | Count-scaled moments in N m |
| `zeroed_fx_N`…`zeroed_tz_Nm` | Corresponding scaled value minus the saved six-element tare vector |

Both count divisors are 1,000,000 in this dataset; each record's actual configuration is saved in `ati_calibration.json`. Tare values, their sample SDs/limits, sample count, observation times, and status are retained. `record_bias` is the vector actually applied to the recording. Sensor axes and mounting directions are retained; these are not anatomical coordinates. The configured six-component tool transform is zero (`cfgtfx`). Positive compression is **minus zeroed Fz**, also stored as `force_down_source_N`; the signed Fz waveform itself is not negated in scalar decomposition. Moments have not been replaced with position channels or transformed into a thumb-contact wrench.

`raw.rdt` retains binary 36-byte ATI packet payloads. `packets.csv` gives byte `offset`, byte `length`, and both receipt timestamps for each payload. `bias_id` indexes the tare, not a participant or specimen. Network addresses, device identities, and private source paths are omitted from extracted calibration metadata.

### WitMotion WT901SDCL IMU

`imu/decoded.csv` and `raw.bin` are byte-identical acquisition copies. The nominal/native sampling rate used by analysis is **200 Hz**. Each timestamped payload occupies 28 bytes; `blocks.csv` records receipt-block byte `offset`, byte `length`, and monotonic/Unix receipt timestamps. A Bluetooth receipt block can contain multiple samples; receipt time is not a separate calibrated sampling time for each frame.

| Columns | Meaning / units |
|---|---|
| `segment`, `raw_offset` | Acquisition segment index and sample byte offset in raw.bin |
| `device_datetime_unmapped` | Device calendar timestamp string before time correspondence; its calendar date is not a verified experiment date |
| `host_receive_monotonic_ns`, `host_receive_unix_ns`, `host_elapsed_sec` | Receipt timestamps and common-gate elapsed time, as for ATI |
| `raw_ax_count`…`raw_az_count` | Signed accelerometer counts; decoded as count / 32768 × 16 g |
| `raw_wx_count`…`raw_wz_count` | Signed angular-velocity counts; decoded as count / 32768 × 2000 deg/s |
| `raw_roll_count`, `raw_pitch_count`, `raw_yaw_count` | Signed Euler-angle counts; decoded as count / 32768 × 180 degrees |
| `ax_g`, `ay_g`, `az_g` | Acceleration in g, in the device axes |
| `wx_deg_s`, `wy_deg_s`, `wz_deg_s` | Signed device-axis angular velocity in deg/s; the paper's wx/wy/wz |
| `roll_deg`, `pitch_deg`, `yaw_deg` | Decoded Euler angles in degrees |

`processed/records/*/imu_orientation.csv` is a derived acquisition-orientation table, not raw measurement. It retains receipt/device fields, `valid`, `relative_qx/qy/qz/qw`, relative Euler angles, `rotation_angle_deg`, and `euler_singular`. The recorded convention is `R_world_imu = Rz(yaw) Ry(pitch) Rx(roll)` and `R_initial_current = R_world_initial.T @ R_world_current`. Quaternions are **XYZW** and relative orientation is expressed in the initial IMU frame. `euler_singular` flags the Euler representation, not a zero-valued motion sample. The orientation reference and measured calibration-quality summaries are provided; no additional anatomical alignment is asserted.

### NOKOV optical exports

Nominal rate: **90 Hz**, confirmed by vendor export headers. `Tracker0` is the fixed reference rigid body; `Tracker1` is the moving wrist-shell rigid body. The measured moving object is the wrist shell, not the thumb joint itself.

`Tracker0.csv` and `Tracker1.csv` retain the complete hierarchical vendor format. The data row containing `Timestamp` marks the start of the actual column header; use `read_pose_csv` in `code/optical_geometry_helpers.py`, not a plain single-header CSV reader.

| Vendor fields | Meaning / units |
|---|---|
| First data column (`Frame#` header spans rows) | Native frame number |
| `Timestamp` | Optical epoch timestamp in milliseconds |
| `XToGlobal1`, `YToGlobal1`, `ZToGlobal1` | Rigid-body position in optical world coordinates, mm |
| `QxToGlobal1`, `QyToGlobal1`, `QzToGlobal1`, `QwToGlobal1` | XYZW orientation, local rigid-body axes to optical world axes |
| `ExToGlobal1`, `EyToGlobal1`, `EzToGlobal1` | Vendor Euler-angle fields; ZYX order is declared in the header; numerical unit not independently verified here and these columns are not consumed |
| `Vx/Vy/Vz/VrToGlobal1`, `Ax/Ay/Az/ArToGlobal1`, `EVx/EVy/EVzToGlobal1`, `EAx/EAy/EAzToGlobal1`, `Segmentlength1` | Unused vendor-derived fields; detailed units/conventions unknown in this release. The analysis calculates its own derivatives. |

TRC files contain frame, native `Time` in seconds, epoch `Timestamp` in milliseconds, and marker XYZ in mm. `Tracker0.trc` and `Tracker1.trc` contain measured Marker1–4 and vendor `RCalc::Marker*` fields. The measured marker fields are distinguished from reconstructed ones; reconstructed coordinates do not replace occluded measured markers. `Unnamed.trc` retains vendor unassigned-marker data, which are not used by the paper. Marker count is declared in each TRC header. Empty coordinates remain missing. Only the private path on the first TRC header line was replaced with the public filename; all subsequent bytes are unchanged. This release provides optical exports, not raw camera images.

## Time coordinates and analysis support

`processed/config/recording_time_origins.csv` retains the common software gate in monotonic and Unix nanoseconds. The original receipt timestamps remain separate from all mapped times. IMU native elapsed seconds are computed relative to the first retained device timestamp. ATI native counter seconds are `(ft_sequence - first_ft_sequence) / 7000` in this cohort; the original fit did not apply counter reset correction. Optical native seconds are the TRC **Time** column, preserved as `optical_s`.

`ati_imu_alignment.csv` contains the frozen device-to-common-software affine coefficients. For IMU, `h = imu_device_to_host_receive_scale * tI + imu_device_to_host_intercept_s`; for ATI, use `ati_counter_to_host_receive_scale` and `ati_counter_to_host_intercept_s`. `optical_clock_models.json` maps native optical time as `h = alpha * tN + beta`. Its arrays give the saved anchor times; residuals are software correspondence diagnostics. They are not hardware latency or physical synchronization accuracy. The release does not refit any map. `leave_anchor_out_max_error_s` is retained with its historical field name; source recovery noted that for two records this value was an in-sample residual, so it must not universally be interpreted as independent held-out error.

`sample_epoch_phase_calibrated`, `physical_sample_latency_verified`, `whole_record_motion_validation_passed`, and `clock_global_motion_certified` remain false. File completeness does not change these statements. Recording/receipt dates and device calendar dates are preserved as measurement metadata; coded identifiers do not remove every identifying feature.

`task_boundaries.csv` contains the adopted `start_optical_s`, `end_optical_s`, their common-host coefficients, task duration and boundary provenance/uncertainty. These saved intervals are the analysis selection, not a new synchronization estimate. `record_index.csv` contains corresponding `start_host_s` and `end_host_s`. All current records have zero saved pause duration; the release uses those frozen intervals instead of rerunning historical stage annotation.

IMU-cycle progress is dimensionless `u = 0, .01, ..., 1` (101 points). For a saved interval `[start_device_s, end_device_s]`, queries are `start + u * duration`; corresponding force/moment queries use the frozen software time correspondence. Native unfiltered physical-channel values define the waveforms. Valid interpolation requires valid neighboring samples no more than **20 ms** apart, and accepted cycles require at least **90%** progress coverage per modality. Unsupported bins stay missing. Recording means weight accepted cycles equally at each supported progress point; `n_cycles` is the contributing count and `q10/q90` describe across-cycle values. No zero padding is introduced for missing support.

`formal_ati_imu_cycles.csv` retains native/host endpoints, durations/frequencies, modality inclusion flags and coverage fractions. `frozen_event_thirds.csv` assigns each cycle by its **host-time midpoint** to one of three equal-duration parts of the adopted formal interval. `third=1,2,3` means early/middle/late. `formal_third_imu_profiles.csv` stores those group means and support counts. Only angular-velocity third-specific waves are included in the manuscript's temporal decomposition.

The IMU detector's saved configuration specifies fourth-order zero-phase Butterworth **1–5 Hz**, 200 Hz, used for event detection, not replacement of native waveform values. The fixed signed axis, fit/validation record IDs, ranges and weights are in `imu_detector.json` and `imu_axis_fit_runs.csv`. Fitting used host time 15 s to each saved endpoint, continuous valid runs with at least 600 samples, uncentered second moments `Z.T @ Z / n`, equal run weights within each record and equal weights for the two fitting records. The largest absolute axis coefficient is positive. Prominence, same-phase spacing and accepted native period limits remain those saved in the configuration.

## Processed native arrays and optical coordinates

NPZ files use named NumPy arrays and require no pickle. File-specific shapes and dtypes are in the schema inventory. Let Nf, Ni and No denote ATI, IMU and optical sample counts for a recording; they need not match.

- `ati_imu_native.npz`: `force_source_row` and `imu_source_row` are zero-based decoded CSV row indices. `force_sensor_zeroed` is Nf × 6 in Fx,Fy,Fz,Tx,Ty,Tz order (N,N m); `imu_gyro_deg_s` is Ni × 3 in device XYZ (deg/s); `imu_relative_quaternion_xyzw` is Ni × 4. Other time/count/validity arrays are one value per corresponding sample. Receipt and device elapsed times are separate. `imu_batch_id` identifies a shared receipt block. `force_down_source_N = -force_sensor_zeroed[:,2]`.
- `force_valid` requires acquisition valid=1 and six finite zeroed values. `imu_valid` requires source orientation validity, finite relative quaternion and finite three-axis gyro. Source orientation validity marks the first sample valid and then device increments in (0,15.001 ms]. `*_source_valid` retains that earlier flag. `force_edge_contiguous` additionally requires two valid samples, RDT increment 1 and increasing host time; `imu_edge_contiguous` requires two valid samples and device increment in (0,5.001 ms]. Their first element is false. These are native continuity diagnostics; frozen phase-level inclusion is stored separately.
- `ati_time_axis.npz`: `ati_ft_sequence`, `ati_device_counter_elapsed_s`, `ati_host_receive_elapsed_s`, `ati_valid`, and `ati_source_row` retain the original counter timeline and row correspondence.
- `optical_geometry_fixed/moving.npz`: N × 3 world positions, N × 4 XYZW rotations, N × 4 marker residuals, per-frame observation/time/geometry flags, a 4 × 3 local marker template in mm, and the actual `template_train_frames`. The template was the componentwise median of the first 180 qualifying four-marker poses within the first 10 native seconds. These geometry inputs are not themselves the final formal-task mask.
- `optical_qualified_poses.csv.gz`: native frame, `optical_s`, final archived `valid`, relative XYZ, relative quaternion, fit residuals, observed-marker counts and diagnostics. The relative transform is `R01 = RW0.T @ RW1`, `p01 = RW0.T @ (pW1 - pW0)`. Position stays in Tracker0 axes. Quaternions represent moving-shell axes relative to Tracker0.
- `optical_native_derivatives.npz`: `task_mask`, `formal_valid_mask`, `finite_spatial_omega_mask`, native frame/time, N × 3 linear velocity (mm/s) and spatial angular velocity (deg/s). Derivatives are computed within valid frame-consecutive runs, never across gaps. Spatial angular velocity is expressed in Tracker0 axes.

Optical diagnostics preserve the original reason bits: **1** invalid exported pose, **2** invalid time continuity, **4** large orientation step, **8** insufficient/nondegenerate marker support, **16** marker geometry residual failure, **32** annotated exclusion. The final relative-pose validity uses the archived quality criteria and the independent rigid fit: at least three noncollinear visible markers, triangle-sine threshold 0.1, rigid-fit RMS at most 3 mm and maximum visible edge-length error at most 4 mm. The jump and label-swap-suspicion fields remain diagnostics; their presence alone is not proof of marker identity failure or a direction to recompute a stricter mask. The packaged workflow consumes the final archived `valid` column.

Optical PC1 is fitted to centered signed spatial angular velocity in each formal record with equal valid-sample weights (sample covariance). Its direction has positive Tracker0 Y. Detection uses native angular velocity, no added 1–5 Hz IMU filter: peak distance 20 frames, prominence `max(.35*SD, 1 deg/s)`, period 0.2–0.7 s. All cycle samples must be valid, formal and consecutive; maximum gap is 1.5/90 s. `formal_native_cycle_candidates.csv`, `formal_cycle_QC.csv`, and axis diagnostics preserve events, inclusion status and fitting output.

For a spatial cycle, translation is relative to its own first position and rotation to its own first orientation. Position interpolation and quaternion SLERP stay inside that valid native cycle. Recording means average cycle positions arithmetically and rotations by the proper chordal rotation mean. `q` and `-q` are equivalent rotations. Logarithmic rotation vectors use the principal SO(3) logarithm and cycle-start moving-frame axes; they are not arithmetic Euler differences.

`optical_cycles/{record_id}.npz` has `array` shape (accepted cycles,101,9): dx/dy/dz mm, vx/vy/vz mm/s, spatial wx/wy/wz deg/s; `quat` shape (cycles,101,4) XYZW, `meanq` (101,4), and `times` (cycles,3) = native start/end/duration seconds. `task_start`/`task_end` are scalars. These arrays were regenerated with the original optical functions, with recording means checked against frozen outputs. `task_profiles` stores `p` (101,3) mm, `q` (101,4), `w` (101,3) deg/s and `n` (101) observed counts. Each whole-task bin requires at least five samples. Whole-task position and orientation use the first formal second's mean position and rotation; this differs from the cycle-start reference. Distances depend on these references.

## Quantitative result fields

For each scalar pair, finite common progress points are equally weighted. `mean_a/b` are centers η; `amplitude_a/b` are centered RMS amplitudes a. Standardized signed shapes are `(x-η)/a`. `distance_squared` is the mean squared waveform difference and exactly decomposes into:

```text
center_MSE    = (ηa − ηb)^2
amplitude_MSE = (aa − ab)^2
shape_MSE     = aa * ab * Ds^2
```

`distance` is its square root, in the channel's physical unit. `shape_RMS_dimensionless` is Ds; `rho` is the mean standardized-shape product. `*_fraction` divides each component by distance squared; fractions are dimensionless contributions, not accuracies or percentages of loading. `unit` and `squared_unit` distinguish N/N², N m/(N m)² and deg/s/(deg/s)². `delta_*_b_minus_a` preserves the explicit direction. `identity_error` is the numerical sum residual. `n_common`, `n_union`, `n_shared_nonfinite`, and `phase_ids` describe exact support. `status` records any undefined quantity; empty fractions/shapes are not replaced by zero. Near-zero handling uses the original 64 × machine-epsilon × magnitude threshold in `decompose`.

`decomposition_summary.csv` reports pairwise-statistic medians, minima and maxima. `largest_*_count` requires a **unique strict maximum** among the three MSE terms. Fig. 4a separately takes the arithmetic mean of per-pair fractions. `record_features.csv` uses each recording's own finite support for its center/amplitude, not another recording's common pair support; optical ranges are native-frame P95−P05 and orientation excursion is the P95 geodesic angle from the specified reference.

`material_cell_means.csv` averages two recording features arithmetically. `material_adjacent_changes.csv` follows P1→P3→P2→P4. `material_monotonicity.csv` records the three adjacent signs, tolerance and increasing/decreasing/mixed/flat category for repetitions and their mean. `material_endpoint_changes.csv` is P4−P1, which does not establish monotonicity of the full sequence. All rows remain stratified by practitioner and round.

Temporal rows use middle−early, late−middle, or late−early contrasts. Each contrast has 32 records × 3 channels. Same-third repeat rows contain 16 pairs × 3 thirds × 3 channels; the reference set is **48 pairs per channel**, not 32 records. Table I reports medians for Δη, Δa and Ds; +/− are counts of positive/negative amplitude changes. `amplitude_*_count` counts records, not samples. Amplitude changes have deg/s units; standardized-shape distances are dimensionless and their magnitudes are not compared across units.

Optical pair tables contain RMS Euclidean path distances in mm and RMS SO(3) geodesic angles in degrees, with axiswise RMS values where supplied. `cycle_mean_position_shape_RMS_mm` is a distance between mean cycle paths, not a spatial-center displacement. Task pairs use common supported task bins (at least 50), distinct from cycle-progress support. Published 19.37/4.34 mm values are medians of the practitioner/repetition pair distances.

## Independent material data

The public raw material CSVs are explicitly **column subsets** of the original recorder files, with numeric text and row order retained. `replay_measured_columns.csv` contains Unix `timestamp` seconds, replay `state`, formal `force_time` seconds and positive-compression `controller_measured_z` in N. `material_measured_columns.csv` contains Unix time, `replay_state`, `fit_input_valid` and `d_rel_m` in metres. The plot uses TRACKING rows and valid material input, the original measured-force interpolation helper and a 20 ms maximum support gap; its display grid is `arange(.25,57.9,.04)` seconds.

`d_rel_m` is relative indentation along the saved surface normal from the shared **stable approximately 2 N preload** reference, not penetration from an unloaded surface. Reference JSON retains base-frame position/sphere-center metres, XYZW quaternion, sphere offset metres, normal vector, measured N, event times and stability criteria. In its `criteria`, target/tolerance/span are N; hold/max_gap/timeout are seconds; position_span is metres; speed is m/s. `max_depth_mm` is a saved configured protection value, not a measured depth.

`results/material/trial_slopes.csv` retains each saved 12–18 N stiffness estimate in N/mm. Its medians reproduce the specimen order; the fitting of those slopes was not rerun. Fig. 1c curves were independently regenerated from measured recorder columns and compared with the frozen display curves.

## Missing values and unverified scope

CSV empty/NaN fields and NumPy NaNs denote unavailable numeric values; boolean/0–1 validity arrays determine eligibility. Zero is a physical numeric value unless an explicit field definition says otherwise. Raw invalid samples remain present; no forward fill, surrogate force channel or new exclusion was introduced.

The runnable starting point, missing standalone ATI/IMU individual-cycle exports, unrerun raw optical qualification, and unrefitted clocks/slopes are stated in the README and verification report. Unknown unused vendor Euler/derivative units remain unknown. No data completeness or numerical agreement claim certifies absolute phase accuracy, anatomical calibration, statistical equivalence, or clinical effects.
