# Discrete-time (slotted) models

Models on a slot lattice rather than on the continuous time axis. Time
advances one slot at a time; in each slot a job in service completes with
probability p and an arrival occurs with probability b, both recorded at the
end of the slot with the departure resolved before the arrival (Daduna's LA
rule and D/A rule). Every rate is therefore a per-slot probability and every
time is a number of slots.

`SolverNC` enters this route only when `options.config.slotted` is true, the
same switch `SolverLDES` uses; it is never inferred from the presence of a
`Geometric` distribution. When the switch is on and the model falls outside
the discrete-time product form, the solver errors instead of approximating.

| Example | Feature | Reference |
|---|---|---|
| `dt_geogeo1` | Geo/Geo/1, unbounded buffer | Daduna (2001), Thm 2.3 / Cor 2.7 |
| `dt_geogeo1_loss` | Geo/Geo/1/L loss system, blocking probability | Cor 2.8 |
| `dt_bernoulli_loaddep` | load-dependent Bernoulli server, arrival theorem | Ex 2.10, Thm 2.11 |
| `dt_cycle` | closed cycle of Bernoulli servers | Cor 3.4, Prop 3.18-3.19, Cor 3.20 |
| `dt_cycle_loaddep` | closed cycle with state dependent service | Thm 3.2 |
| `dt_cycle_multiclass` | multichain closed cycle | Sec 3.2 |

Reference: H. Daduna, *Queueing Networks with Discrete Time Scale*, LNCS 2046,
Springer, 2001.
