# Progress Log

## 2026-07-16

- Set up repo structure: `src/`, `include/`, `data/`, `figures/`,
  `analysis/`, `docs/`, `report/`, root `CMakeLists.txt`.
- Initialized git repository, initial commit.
- No tasks started yet. Next: Task 1 (definitions/notation, no code) as
  a warm-up, then Task 2 (π estimation).

## 2026-07-17

- Task 2 (π estimation) implemented end-to-end: `src/pi_estimator.cpp`
  runs both the area method and the MC-integration method, writing
  convergence data (`data/task2/pi_estimation_convergence.csv`) and a
  fixed-N distribution over M experiments
  (`data/task2/pi_estimation_distribution.csv`).
- `src/pi_estimator.cpp` has **not been code-reviewed yet** — still to
  go through it for idiomatic C++ before considering Task 2 done.
- Set up the report/figures pipeline: plotting scripts in `figures/`
  now `savefig` PDFs into `report/figures/<task>/`, and
  `\graphicspath{{figures/}}` was added to `template.tex`/`task2.tex`
  so `\includegraphics` can use short relative paths. Verified with a
  clean `latexmk -pdf task2.tex` build.
- `figures/task2.py` (convergence plot) and `figures/task2_3.py`
  (fixed-N histogram, area vs. integration) both run and save
  correctly now, after fixing two bugs along the way:
  - an f-string with clashing quote characters
    (`f"...{df["method"][0]}..."`) — a SyntaxError on Python < 3.12
    (PEP 701 relaxed this in 3.12+, but worth not relying on it);
  - `df_to_plot["method"][0]` used label-based indexing into a
    boolean-masked slice, which fails with `KeyError` whenever the
    slice doesn't happen to contain a row labeled `0` — fixed with
    `.iloc[0]` for positional access.
- `report/task2.tex` currently has the convergence figure wired in
  with a real caption; the distribution figure is present but only a
  bare `\includegraphics`, no discussion text yet — report is still
  barebones (plots dropped in, no written analysis/interpretation).

### Next / open for Task 2

- Write up the actual analysis text in `task2.tex`: describe both
  methods, compare mean/variance between area vs. integration (numbers
  already printed by `task2_3.py`, just need to be pulled into prose
  or a table), and comment on the convergence behavior shown in the
  plot.
- Finish the distribution figure's caption/label and surrounding
  discussion.
- Go back through `src/pi_estimator.cpp` for a proper code review
  (idiomatic C++ pass) before calling Task 2 done.
- Then move on to Task 3 (CDF inversion sampling).

## 2026-07-18

- Task 3 (CDF inversion) implemented in `src/cdf_inversion.cpp`:
  samples τ via CDF inversion, writes samples to
  `data/task3/cdf_inverse_data.csv` for the histogram, and computes
  the MC estimates of I1/I2 with error bars from sample variance.
  Went through a code review pass (idiomatic C++: pass-by-value
  callable instead of a forwarding reference, factored out
  `sample_tau` to avoid duplicating the CDF inversion formula, still
  need to hoist the normalization constant `1-exp(-5)` out of its
  remaining duplicated spots).
- `figures/task3_1.py` plots the histogram (normalized via
  `density=True`) against the exact curve on a single shared y-axis.
- Added `data/*/*.csv` and `report/figures/` to `.gitignore` (data is
  regenerated non-deterministically each run, no fixed seed) and
  untracked the previously-committed `data/task2/pi_estimation_convergence.csv`.
- Report for Task 3 not started yet.

### 2026-08-09

- Task 4 was implemented to completion
- Did a blocking analysis method and optimized step size delta to minimize the autocorrelation time
- Results match and variance from the blocking analysis looks correct
- Report finished in its raw form (some improvements to be done later)

### 2026-08-multiday

- Task 5 is fully implemented
- Task 6 is also nearly completely implemented, but is missing some fine tuning and improvement which is the current goal.