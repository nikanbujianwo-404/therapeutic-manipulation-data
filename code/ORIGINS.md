# Code provenance and release adaptation

| File | Original functional source | Release changes |
|---|---|---|
| `analyze_quantitative.py` | `analyze_quantitative_findings.py`: exact `decompose` function and current-paper feature, pair, material and temporal loops | Relative output root; public record labels; unused historical analyses removed. Numerical rules retained. |
| `optical_core.py` | `reanalyse.py`: `runs`, `meanR`, `derivatives`, `profile`, `cycle_profiles`, `formal_native_events`, `compare_pairs` | Exact function extraction with the original required constants/imports. |
| `optical_geometry_helpers.py` | `reconstruct_spatial.py`: vendor CSV reader and rigid-fit helpers | Exact function extraction; not called by the reproduction entry and not a complete raw geometry pipeline. |
| `phantom_force_signal.py` | Original measured-force interpolation helper of the indentation experiment | Byte-identical copy. |
| `plot_figures.py` | Current scientific plotting lineage, with manuscript-approved public expert labels | Relative paths; public record IDs; numerical Fig. 1c only; Matplotlib 3.5 `ncol` keyword. Data selections and numerical operations retained. Later publication panel cropping is not reproduced. |
| `reproduce.py` | New portable orchestration around those kernels and frozen current inputs | No new fitting/clock estimation; writes a new output directory; compares rows by scientific identity because public IDs sort differently. |
| `validate_data.py` | New release integrity/readability validation | Checks manifest design, all file hashes, readable arrays, raw-to-native sample values, clocks and recording coverage. |

Numerical output comparisons, rather than filename ordering, establish agreement with the frozen results. Source hashes and internal source locations are retained in the private handover, not in public metadata. Dependency versions are those actually used in release validation. The wrapper invokes child Python processes with `-s` to prevent unrelated user-site packages from overriding the tested environment.
