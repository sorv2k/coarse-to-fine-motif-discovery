# Coarse-to-Fine Search for Temporally Scaled Time Series Motif Discovery

MS Capstone Project, UC Riverside, March 2026  
Author: Sourav Guruprasad  
Advisor: Professor Eamonn Keogh

## Overview

Given two time series A and B, this project finds a subsequence in A 
and a subsequence in B that match when one is rescaled in time by a 
factor s. This is useful when the same pattern appears at different 
speeds, for example ECG signals from a young and elderly patient, or 
insect calls recorded at different temperatures.

The brute force approach tests all 101 scaling factors from 0.5 to 1.5 
at 0.01 steps exhaustively. The Coarse-to-Fine Search tests only 11 
coarse factors at 0.1 steps, identifies the best region, then refines 
with approximately 31 fine tests. Total approximately 42 tests vs 101 
for brute force.

## Results

| Dataset | Brute Force | Coarse-to-Fine | Speedup |
|---|---|---|---|
| 20,000 points (synthetic) | 171 sec | 70.87 sec | 2.8x |
| 200,000 points (synthetic) | 286 min | 31 min | 9.2x |
| ECG Fantasia (real data) | 0.779 expected | 0.790 found | 1.41% error |

Success rate on synthetic validation: 83% (5 of 6 test cases)

## Requirements

- MATLAB R2020a or later
- No additional toolboxes required

## Files

| File | Description |
|---|---|
| `mpx.m` | Matrix Profile function, core dependency |
| `final_FindRescaledMotif_Adaptive.m` | Main Coarse-to-Fine Search algorithm |
| `final_generate_signals_with_scale.m` | Synthetic signal generator with embedded motif |
| `final_test_robustness.m` | 6-case synthetic validation with plots |
| `final_smoothness.m` | Distance curve smoothness analysis across 6 test cases |
| `final_test_long_signal.m` | 200,000 point brute force vs coarse-to-fine comparison |
| `explore_fantasia_ecg.m` | ECG extraction and validation on Fantasia Database |

## Usage

### Run synthetic validation
```matlab
final_test_robustness.m
```

### Run large dataset test
```matlab
final_test_long_signal.m
```

### Run ECG real data validation

Download the Fantasia Database from PhysioNet and update the 
csv_path variable in explore_fantasia_ecg.m, then run:
```matlab
explore_fantasia_ecg.m
```

### Use the algorithm directly
```matlab
result = final_FindRescaledMotif_Adaptive(A, B, m, 0.5, 1.5, 0.1, false);
fprintf('Best scaling factor: %.3f\n', result.bestScale);
```

Parameters:
- A: time series A (vector)
- B: time series B (vector)
- m: motif length
- MinSF: minimum scaling factor (e.g. 0.5)
- MaxSF: maximum scaling factor (e.g. 1.5)
- StepS: step size, use 0.1 for coarse-to-fine, 0.01 for brute force
- PlotFlag: true to show plots, false otherwise

## Algorithm

Phase 1 (Coarse Search): Test 11 scaling factors at 0.1 step size 
from 0.5 to 1.5. Identify the approximate region of the minimum.

Phase 2 (Fine Search): Refine within best region plus or minus 0.15 
at 0.01 steps, approximately 31 tests.

Total: approximately 42 tests vs 101 for brute force.

## Dataset

Real data validation uses the Fantasia Database from PhysioNet:
- 2-hour ECG recordings at 250 Hz
- 20 young subjects (21-34 years) and 20 elderly subjects (68-85 years)
- Citation: Iyengar N, et al. Age-related alterations in cardiac 
  interbeat interval dynamics. Am J Physiol 1996.

## References

- Yankov D, Keogh E, et al. Detecting time series motifs under uniform 
  scaling. KDD 2007.
- Yeh CM, Keogh E, et al. Matrix Profile I: All pairs similarity joins 
  for time series. IEEE ICDM 2016.
