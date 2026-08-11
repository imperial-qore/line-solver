# libQBD - MATLAB port

MATLAB port of [libQBD](https://github.com/ProgGrey/libQBD), a library for
analyzing Quasi-Birth-Death (QBD) processes.

**Original author:** Sergey Astaf'ev, IAMR Karelian Research Centre RAS
**License:** BSD 3-Clause (see LICENSE)

## Features

- Stationary distribution via logarithmic reduction
- Transient analysis of state probability distributions via Taylor series (classic and adaptive)
- Mean number of clients and queue length computation

## Usage

```matlab
addpath(fileparts(mfilename('fullpath')));

proc = libqbd.QBD();
proc.add_zero_level(lambda);
proc.add_level(mu, lambda);
proc.add_final_level(mu);

model = libqbd.StationaryDistribution(proc);
rho = model.get_rho();
L = model.get_mean_clients();
```
