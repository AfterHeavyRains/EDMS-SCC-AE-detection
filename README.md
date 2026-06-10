# Entropy-Driven Hierarchical AE Detection for SCC

Code and supporting outputs for the paper *Entropy-Driven Hierarchical Framework for Robust Acoustic Emission Detection with False Alarm Control*.

The MATLAB scripts produce the entropy/time–frequency figures and the evaluation, ablation, and SNR results; the Python scripts produce the summary and trade-off plots. Generated figures and table summaries are kept under `results/` so the outputs can be checked without rerunning everything.

## Layout

```
code/
  matlab/figure_scripts/   representative signal / entropy figures
  matlab/evaluation/       parameter sweep, evaluation, ablation, SNR
  python/plotting/         summary plots and trade-off maps
data/                      a couple of example samples per thickness + phase map data
results/                   generated figures, tables, and SNR outputs
requirements.txt
```

## Scripts behind the main results

| Item | Script |
|---|---|
| Fig. 2 (CWT/SWT contour) | `code/matlab/figure_scripts/Tensile_250211numericial.m` |
| Fig. 4 (Tsallis entropy) | `code/matlab/figure_scripts/Tensile_250214_SWT_Tsallis_lowE.m` |
| Fig. 5 (entropy rate) | `code/matlab/figure_scripts/Tensile_250221_SWT_TF_lowE_rate.m` |
| Fig. 6 (windowing) | `code/matlab/figure_scripts/Tensile_250228tf_three_lowE_window.m` |
| Phase map | `code/matlab/evaluation/Phase_Map_Visualization.m` |
| Detection pipeline + result figures (AGR/EOR, event categories, weak-event quality, Δt) | `code/matlab/evaluation/Entropy_Detection_Viz_updated_all.m` |
| Main evaluation table | `code/matlab/evaluation/eval_and_save_all_basepara.m` |
| Parameter optimization | `code/matlab/evaluation/optimize_and_sweep_*.m` |
| Ablation study | `code/matlab/evaluation/run_ablation_studyBEIJING.m` |
| SNR characterization | `code/matlab/evaluation/compute_SNR_noise_characterization.m` |
| Summary / trade-off plots | `code/python/plotting/*.py` |

## Running

MATLAB (with Signal Processing and Wavelet toolboxes): open the repo root and run
`addpath(genpath(pwd));`, then run the script of interest. Some scripts contain
machine-specific paths near the top that should be pointed at your local copy.

Python: `pip install -r requirements.txt`, then run the scripts under
`code/python/plotting/`.

Figures 2 and 4–6, the phase map, and the Python summary plots run directly
from the example data and the CSV outputs included here. The four-thickness
evaluation, ablation, and SNR scripts need the full waveform dataset (see below);
their generated outputs are already provided under `results/` for reference.

## Data

`data/` contains a few example per-thickness samples and the phase-map `.mat`
used by the figures above. The complete four-thickness waveform dataset is large
and is available from the corresponding author upon reasonable request.
