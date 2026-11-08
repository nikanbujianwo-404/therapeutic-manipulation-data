# Physical Variation in Expert Therapeutic Manipulation Across Materials and Time

Data and analysis code accompanying the manuscript by **Jixian Zheng, Chengdi Hou, Feng Hao, Jun Hu, Pengcheng Wang, and Lifeng Zhu**. Corresponding author: Lifeng Zhu.

This study describes variation in sustained therapeutic manipulation across practitioners, materials, and time. It separates differences between recording-mean cycle waveforms into mean-level, RMS-amplitude, and standardized signed-shape contributions, alongside independently analyzed wrist-shell spatial paths. Same-condition repetitions provide an empirical reference for difference magnitudes.

This repository provides the measurement data, processed analysis inputs, and analysis code accompanying this study. Download the data packages from the release linked below and follow the instructions in “Download and extraction.”

- Release version: `1.0.0`
- Repository: https://github.com/nikanbujianwo-404/therapeutic-manipulation-data
- Data and code reuse permissions: **awaiting author decisions**, separately; see [LICENSE_STATUS.md](LICENSE_STATUS.md).

## Scope and identifiers

The main experiment comprises **32 recordings: two practitioners × four specimens × two rounds × two repetitions**. Public record IDs use `ExpertA` and `ExpertB`, consistently with the manuscript. These are coded identifiers, not a claim of strict anonymization. Four specimens are ordered by their independent mechanical responses as **P1, P3, P2, P4**. Twelve independent indentation trials (three per specimen) are indexed separately and are not counted as manipulation recordings.

The release contains full six-axis ATI measurements, IMU records, NOKOV rigid-body/marker exports, original sample times and quality information, frozen task intervals/events/mean waveforms, optical calibration parameters, and recording/pair-level result tables. Historical cohorts and other research projects are excluded. NOKOV camera images are not included. Apparatus photographs in Fig. 1a,b are not part of this numerical release; the confirmed robotic indentation panel remains in the manuscript.

## Download and extraction

Open the repository's **Releases** page, select `v1.0.0`, and download both data assets. GitHub's automatically generated “Source code (zip)” does **not** replace the data assets.

| Asset | Contents | Extract into |
|---|---|---|
| `TBME_raw_measurements_v1.0.0.zip` | `data/raw/`: 32 multimodal recordings and selected measured columns from 12 material trials | Repository root |
| `TBME_analysis_data_v1.0.0.zip` | `data/processed/` and `data/results/`: analysis inputs, masks, events, calibration/time parameters, cycle paths and result tables | Same repository root |

Each ZIP is independent and uses ordinary ZIP compression. Download `SHA256SUMS.txt` and `FILE_MANIFEST.csv` from the release to verify the assets and extracted data. No external symlinks are required.

```bash
# Run in a new directory containing the repository files and downloaded assets.
sha256sum -c SHA256SUMS.txt
unzip TBME_raw_measurements_v1.0.0.zip
unzip TBME_analysis_data_v1.0.0.zip
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -s code/validate_data.py --data data
python -s code/reproduce.py --data data --output reproduced
```

`SHA256SUMS.txt` also lists `GITHUB_REPOSITORY_FILES.zip`; download that asset as well for the full checksum command, or check the selected ZIP entries separately. Tested with Python **3.10.12** and the dependency versions in `requirements.txt`. Installation on a fresh network-connected system was not tested; the complete workflow was run with those installed versions in an isolated extraction directory. Choose a new output directory for each run. Scripts do not overwrite inputs or an existing reproduction directory.

## Reproduction boundary

The runnable workflow starts from **frozen recording-mean ATI/IMU waveforms and qualified native optical poses**. It recomputes scalar decomposition, material responses, temporal comparisons, optical derivatives/events/cycle paths/task paths and pair distances. It also regenerates the 12 material force–depth curves from measured recorder columns. Numerical panels for Figs. 1c and 2–6 and Table I values are generated; publication-specific photographic composition and later panel cropping are not reproduced.

Raw acquisition files are provided for provenance and further work. The packaged entry does not claim a complete raw-to-paper pipeline: it does not regenerate primary optical marker-quality masks, ATI/IMU events or recording-mean waves from raw samples, estimate clocks, or refit material slopes. Frozen ATI/IMU cycle boundaries, native sample arrays, masks and time parameters are supplied, but per-cycle ATI/IMU query-support arrays were not recovered as a standalone export. Individually reconstructed optical cycle arrays are included and regenerated.

Saved software time correspondence is **not** verified physical synchronization. The original `physical_sample_latency_verified=false`, `whole_record_motion_validation_passed=false`, and `sample_epoch_phase_calibrated=false` statuses are retained. Native time, mapped software time, and dimensionless cycle progress are distinct quantities; see [data_dictionary.md](data_dictionary.md).

## Structure and outputs

```text
record_manifest.csv            32-record index and relative data paths
material_manifest.csv          12 independent indentation trials
data_dictionary.md             definitions, units, time and coordinate conventions
schemas/                      file-level dimensions and per-field descriptions
code/                         portable current-analysis scripts and original kernels
config/                       release metadata
data/raw/records/              original ATI/IMU payloads and decoded tables; NOKOV exports
data/raw/material/             measured recorder column subsets
data/processed/records/        native arrays, masks, relative poses, calibration
data/processed/config/         frozen task boundaries, detector and clock parameters
data/processed/evidence/       recording-mean waves, cycles and temporal means
data/processed/optical/        frozen optical events and quality tables
data/processed/optical_cycles/ individually reconstructed optical cycle arrays
data/processed/task_profiles/  whole-task profiles on 101 bins
data/results/                  frozen current-paper results and figure-source tables
reproduced/                    created by code/reproduce.py
```

| Manuscript evidence | Input / frozen result | Regenerated output under `reproduced/` |
|---|---|---|
| Fig. 1c; specimen mechanical order | raw material recorder subsets; `results/material/trial_slopes.csv` | `material/mechanical_curves.csv`; curves recomputed, saved slopes retained |
| Fig. 2 | `processed/figure_sources/P1_*`; `processed/evidence/formal_ati_imu_phase_profiles.csv` | `figures/Fig2_states_signed_cycles.*`, matching source tables |
| Fig. 3; spatial cycle/task comparisons | qualified native poses and frozen task intervals | `optical/cycle_pair_distances.csv`, `optical/task_pair_distances.csv`, `optical/cycles/`, `optical/task_profiles/` |
| Fig. 4; Fz/wx/wy contributions | recording-mean waves and matched pair index | `evidence/quantitative/pair_decomposition.csv`, `decomposition_summary.csv` |
| Fig. 5; material responses | per-record waveforms and native optical features | `evidence/quantitative/record_features.csv`, `material_*` |
| Fig. 6; early-to-late changes | frozen third-specific IMU mean waveforms | `evidence/quantitative/temporal_*.csv` |
| Table I | same third-specific means | `figure_sources/Fig6_all_contrast_summary.csv` |

`verification/REPRODUCTION.json` records numerical comparisons against all 11 frozen quantitative tables, optical means/events/pair distances, material curves, and the manuscript's highlighted results. Scalar fractions are calculated per pair before their medians; Fig. 4a uses arithmetic mean fractions. Pair counts and cycles are not independent participant counts.

## Citation and permissions

Until a publication DOI exists, cite the manuscript title and authors above and identify the data/code archive by its real version and DOI once assigned. Do not use the placeholders as identifiers. No manuscript PDF, reference-library PDF, third-party library, environment, model checkpoint, or submission correspondence is distributed here.

The code contains extracted original scientific kernels and a new path/verification wrapper. `code/ORIGINS.md` identifies the functional provenance and packaging changes. Third-party dependencies are installed separately under their own licenses. The authors must choose and document appropriate data and code licenses before public publication.
