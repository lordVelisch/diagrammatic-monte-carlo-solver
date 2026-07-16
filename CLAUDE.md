# Project: Basic Diagrammatic Monte Carlo

C++/Python project implementing a sequence of Monte Carlo methods leading
to a basic DMC code, for a computational condensed matter course. Also a
deliberate C++ fluency exercise ahead of a lattice QCD / Kokkos internship
— treat code quality and idiomatic C++ as first-class, not incidental.

Full task list and equations: @docs/tasks.md (task PDF also in `docs/`).
Session-to-session history: @docs/progress.md — read it at the start of a
session before assuming this is a cold start.

## Architecture

- **C++**: all simulation code (RNG, MC/MCMC/DMC loops, updates,
  measurement accumulation). Writes plain CSV to `data/<task>/`.
- **Python** (`analysis/`, numpy + matplotlib): all plotting and
  comparison to analytic curves, reading that CSV. No plotting logic in
  C++.
- Repo layout: `src/`, `include/`, `data/`, `figures/`, `analysis/`,
  `docs/`, `report/` (LaTeX), root `CMakeLists.txt`.

## Build

CMake, out-of-source `build/`. C++17 minimum. Prefer standard library
(`<random>`, `<vector>`, etc.) over external dependencies — the point is
to learn the language.

## Learning mode vs. glue mode — the most important rule

Default to **guiding**, not solving, for core algorithmic/physics work:
the MC loop structure, CDF inversion, Metropolis-Hastings acceptance
logic, the blocking method, the diagram/update abstraction (Task 5), the
add-β/remove-β updates (Task 6), the Task 7 estimator. Ask leading
questions, review and critique code I've written, escalate hints
gradually — don't write full solutions for these.

For glue — CMake, file I/O boilerplate, plotting scripts after the first
of a given type, repetitive LaTeX — just do it efficiently, no back-and-
forth needed.

If a task doesn't clearly fall into either bucket, ask before starting.

## Conventions

- Idiomatic modern C++: RAII, no raw owning pointers without reason,
  `std::vector`/`std::array` over C arrays, `const`/references used
  properly. Explain *why* when suggesting a more idiomatic pattern.
- Commit per logical unit of work with real messages, not "wip".
- Keep the Task 5 diagram/update framework general enough that Task 6 is
  additive, not a rewrite. Flag it if a design choice threatens that.
- Update `docs/progress.md` after each session: what was done, what's
  next, open questions.
