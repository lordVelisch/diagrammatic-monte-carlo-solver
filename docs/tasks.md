# Task Reference

Full equations and details are in `project_basic_dmc.pdf` in this
directory — this file is a working summary for quick lookup, not a
replacement for it.

1. **Definitions/notation** — no code. Probability, PDF/CDF, expectation,
   variance, etc. Skip implementation, just make sure the concepts are
   solid before Task 2.

2. **π estimation** — area method (count points in unit circle) vs. MC
   integration (∫√(1-x²)dx via importance sampling with uniform p(x)).
   Run M experiments at fixed N, histogram the results, compare mean and
   variance of both methods.

3. **CDF inversion** — sample Q(τ,α) = exp(-ατ), α=1, τ∈[0,5], via CDF
   inversion. Histogram vs. exact curve (get normalization right). MC-
   integrate I1 = ∫τQ dτ and I2 = ∫τ²Q dτ over [0,5], with error bars from
   sample variance.

4. **MCMC** — same target as Task 3, now via Metropolis-Hastings with
   p(x) = U(0,5). Track accept/reject counts; discuss how to raise
   acceptance toward 1. Repeat the I1/I2 integrals. Implement the
   blocking method (course notes ref: Gubernatis/Kawashima/Werner Ch.
   3.4) for correct error bars on correlated samples.

5. **DMC I** — reframe Task 4's code as diagram sampling. 0th-order
   diagram = straight line of length τ, weight exp(-ατ). Build a general
   diagram/update representation (this is the key C++ design task).
   Implement "change-τ" update (self-inverse). Then add "change-α",
   selected with equal probability to change-τ; run for α ∈ {0.5, 1.0}.

6. **DMC II** — add the 2nd-order diagram (2 vertices, coupling V=0.5,
   β ∈ {0.25, 0.75}), per Eq. 8/9 in the PDF. Implement add-β (0→2) and
   remove-β (2→0) updates with their acceptance ratios (Eq. 11/12) — pay
   attention to the context factors p_add/p_rem. Compare against the
   analytic Q(τ,α,V). Think about histogram normalization when the exact
   normalization isn't analytically available in general.

7. **(Optional, advanced) Exact estimator** — derive an unbiased MC
   estimator for Q(τ₀,α,V) at a specific τ₀ (removes binning error), per
   Mishchenko et al. Section 3.C. Implement and test against Tasks 5/6.

## Notes to self (add as they come up)

-
