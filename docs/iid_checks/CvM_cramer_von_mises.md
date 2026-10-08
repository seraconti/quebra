# CvM - Cramer-von Mises test for renewal, time-censored

## The equation

The third functional of the same tied-down Brownian bridge that gives C1 and C2. Kvaloy &
Lindqvist (Technometrics 62(1):101-115, 2020, Section 3.3) derive a class of tests by applying
different functionals to the normalised counting process; this one integrates the square
without a weight:

```
integral W0(s)^2 ds            <- CvM
integral W0(s)^2/(s(1-s)) ds   <- C2 (Anderson-Darling)
integral W0(s) ds              <- C1 (Lewis-Robinson)
```

Derived rather than transcribed. The paper prints its eq (6) middle term as
`- i N (T^2_{i+1} - T^2_i)/tau`, which is dimensionally inconsistent with the first term
`i^2 X_{i+1}/tau`. The typo is in the paper itself, in both versions: eq (6), p. 103, of
the Technometrics paper and eq (6), p. 5, of its preprint arXiv:1802.08339v1. The authors' own R code (`CvMtestobs` in
github.com/jtkgithub/trendtests, identical to the journal's supplementary material) divides
by `tau^2`. Integrating from the definition gives `/tau^2`:

```
integral_0^1 (N(s tau) - sN)^2 ds
    = sum_{i=0}^{N-1} [ i^2 (u_{i+1} - u_i) - i N (u^2_{i+1} - u^2_i) ]
      + N^2 [ u_N^2 - u_N + 1/3 ],        u_i = T_i / tau,  u_0 = 0
```

which reproduces the paper's own tail term exactly. The normalised bridge is
`V0(s) = (N(s tau) - sN) / (gamma_hat sqrt(N))`, giving the `(1/gamma_hat^2)(1/N)` prefactor.

At `gamma_hat = 1` the statistic is exactly the classical one-sample Cramer-von Mises `W^2` on
`u_i = T_i/tau`.

## Where it lives

`analyzers/checks/cvm_cramer_von_mises.py`. `statistic` and `statistic_batch` compute it,
`cvm_limiting_cdf` gives the asymptotic p-value, `run` is the entry point.

`tests/test_checks_cvm.py` pins the transcription from four directions: the `gamma_hat = 1`
identity against `scipy.stats.cramervonmises`, the bracket against the bridge integral
computed by quadrature, the limiting CDF against published critical values, and
`statistic_batch` against `statistic` under the identity permutation.

## Why it is in the battery

The paper reverses its own single-process preference for `m > 1`: its simulation study
(Technometrics Section 6, p. 106; arXiv v1 Section 5) leaves Anderson-Darling out for several
processes "as the Cramer-von Mises test had better level properties in this case". Gapped
records are routine here and `m > 1` is the gap case, so CvM is the goodness-of-fit form the
paper keeps for several processes in its own simulations. The paper's own `m > 1` CvM and AD
statistics are sums weighted by `tau_j` with a simulated null (Section 8.2, Table 4); the
`m > 1` route here is an UNWEIGHTED sum with a permutation null. The paper does not propose
the unweighted sum; it suggests a permutation null only in general terms (Section 9, p. 113).

It also survives a case C2 cannot. On the in-spec clock of a carved record `tau == T_N`, which
makes eq (7) singular and silences C1 and C2; CvM's integrand carries no `1/(s(1-s))` weight,
so its statistic stays finite there.

## Its rows

Two of the nine in `battery.ROW_KEYS`: one asymptotic, one permutation-calibrated. The
asymptotic row runs only when `include_tau_checks` is True and the record has one segment,
as C2's asymptotic row does; C1's needs only the flag. Its limiting null still assumes a
truncation time chosen independently of the events, and an event-determined `tau` breaks
the tied-down bridge for CvM exactly as it does for the other two. Only the permutation row
is entitled to that case.

`ROW_KEYS` is the schema of the bench tables, so registering a check without re-running the
bench makes `bench_acceptance_at_n` return None and every ledger row read
`underpowered / no bench cell`.

## Limitations

Shared with C1 and C2: the asymptotic calibration divides by an estimated `gamma_hat`, so its
finite-n level misses nominal by an amount and in a direction that depend on the gap
distribution. `test_tier3_calibration.py` measures that level and asserts only a wide bound,
deliberately. `LIMITATIONS.md` covers what all six checks share.

The permutation row has no such caveat - it is exact under exchangeability - but it needs a
declared seed, and `run_battery` raises rather than defaulting one.
