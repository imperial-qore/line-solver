# MAMSolver ETAQA

This directory contains the ETAQA (Efficient Truncation and Aggregation for Queueing Analysis)
algorithm implementation from MAMSolver, developed at the College of William & Mary.

## Source

- **Website**: https://www.cs.wm.edu/MAMSolver/
- **Contact**: MAMSolver@cs.wm.edu

## Description

The ETAQA algorithm provides efficient computation of stationary probabilities and
queue length moments for structured Markov chains:

- **QBD** (Quasi-Birth-Death) processes
- **M/G/1 type** Markov chains
- **GI/M/1 type** Markov chains

## Files

| File | Description |
|------|-------------|
| `MG1_G_ETAQA.m` | Computes G matrix for M/G/1 type chains |
| `MG1_pi_ETAQA.m` | Computes aggregated stationary probabilities for M/G/1 type |
| `MG1_qlen_ETAQA.m` | Computes queue length moments for M/G/1 type |
| `GIM1_R_ETAQA.m` | Computes R matrix for GI/M/1 type chains |
| `GIM1_pi_ETAQA.m` | Computes aggregated stationary probabilities for GI/M/1 type |
| `GIM1_qlen_ETAQA.m` | Computes queue length moments for GI/M/1 type |
| `QBD_G_ETAQA.m` | Computes G matrix for QBD processes |
| `QBD_pi_ETAQA.m` | Computes aggregated stationary probabilities for QBD |
| `QBD_qlen_ETAQA.m` | Computes queue length moments for QBD |
| `ParseOptPara.m` | Utility for parsing optional parameters |

## Dependencies

Requires the SMCSolver library by Prof. Benny Van Houdt:
- QBD files: http://win.uantwerpen.be/~vanhoudt/tools/QBDfiles.zip
- MG1 files: http://win.uantwerpen.be/~vanhoudt/tools/MG1files.zip

These are included in LINE under `matlab/lib/thirdparty/QBDfiles/` and `matlab/lib/thirdparty/MG1files/`.

## Usage in LINE

LINE uses MAMSolver ETAQA for solving BMAP/M/1 queues via the `solver_mam_bmap_m_1` function,
which models batch Markovian arrival processes with exponential service as M/G/1 type chains.

## References

When using this code, please cite:

- V. Stathopoulos, A. Riska, Z. Hua, and E. Smirni. "ETAQA Solutions for
  Infinite Markov Processes with Repetitive Structure." *SIGMETRICS Performance
  Evaluation Review*, 39(4):60–72, 2012.

- A. Riska and E. Smirni. "ETAQA: An Efficient Technique for the Analysis of
  QBD-Processes by Aggregation." *Performance Evaluation*, 54(2):151–177, 2003.

## License

This software is provided by the College of William & Mary. Please contact
MAMSolver@cs.wm.edu for licensing inquiries.
