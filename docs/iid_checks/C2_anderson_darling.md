# C2 - Anderson-Darling test for renewal, time-censored

## The equation

Kvaloy & Lindqvist, Technometrics 62(1):101-115 (2020), eq (7), p. 104. The same
equation is eq (7), p. 5, of the Technometrics paper's preprint, arXiv:1802.08339v1.

```
AD = (1/gamma_hat^2) * (1/N) * {
        sum_{i=1}^{N-1} [ (N-i)^2 * ln((tau - T_i)/(tau - T_{i+1}))
                          + i^2 * ln(T_{i+1}/T_i) ]
        + N^2 * [ ln(tau/(tau - T_1)) + ln(tau/T_N) - 1 ]
     }
```

With `u_i = T_i/tau`, eq (7) at `gamma_hat = 1` is exactly the classical Anderson-Darling
statistic `A^2 = N * integral (F_N(u) - u)^2 / (u(1-u)) du`, testing the event times for
uniformity on `[0, tau]`. The `1/gamma_hat^2` generalises it from Poisson to renewal.

## Where it lives

`analyzers/checks/c2_anderson_darling.py`. Two tests pin the transcription:

- `test_eq7_is_the_classical_anderson_darling` compares it against the textbook
  `-N - (1/N) sum (2i-1)[ln u_i + ln(1 - u_{N+1-i})]` form (Anderson and Darling, JASA 49:765-769, 1954, eq (2), p. 765). Two unrelated algebraic forms
  agreeing to float precision is strong evidence the transcription is right.
- `test_eq7_at_gamma_one_matches_the_limiting_ad_null` checks its rejection rate against
  Marsaglia & Marsaglia's `adinf` over 6000 replicates.

**Both force `gamma = 1`, so neither pins the shipped path.** That distinction matters and
is load-bearing - see the limitations below.

## Limitations

**The shipped asymptotic path is oversized at small n, and by more than the transcription
tests can see.** Dividing by an ESTIMATED `gamma_hat` fattens the upper tail: measured at
40,000 replicates on exponential gaps with the shipped sample divisor, 0.0557 at n = 20
and 0.0502 at n = 50 against a nominal 0.05, where the `gamma = 1` path sits at 0.0499.
`test_shipped_c2_asymptotic_size_matches_its_measurement` records those figures and their
seed, and holds the shipped size to a coarse band around them; it is too wide to detect the
oversize itself. The population divisor gives 0.0659 and 0.0533 on the same draws, because
C2 carries `1/gamma_hat^2` and so moves by `N/(N-1)`.

**`tau == T_N` makes it singular.** Its `i = N - 1` summand, `ln((tau - T_{N-1})/(tau - T_N))`, is `+inf` when the last
event lands exactly on the truncation time. That is failure censoring, a different sampling
scheme from the time censoring the equation is derived for, so the check declines it. On the
in-spec clock of a carved record this is the common case, not a corner one: ~73% of
synthetic replicates, and every threshold above 3 µs on the real 0704 record.

Float summation order matters here. `np.sum` is pairwise and `np.cumsum[-1]` is sequential;
they disagree by ulps in either direction on 23% of random vectors. A guard written against
one and a statistic written against the other let a segment pass validation and then produce
`inf`, reported as **p = 0.0** - a false rejection. Both now compare against
`last_event_time` with a documented relative margin.

**There is no asymptotic calibration for m > 1.** The paper's Section 4.2 says the normal
approximation works "less well for the Anderson-Darling test due to the very skew
distribution", and its simulation study (Section 6, p. 106; Section 5 of arXiv v1) leaves AD
out for m > 1 "as the Cramer-von Mises test had better level properties in this case". Asking
for it raises: the paper reports that the approximation works less well and leaves it out of
its m > 1 study, although the authors' R code does compute it. The paper does still
report AD for m = 6 in Section 8.2 (Table 4, p = 0.076), with weights proportional to
`tau_j` and a null simulated from the limiting law, not a normal approximation.

## What we do

Use the permutation calibration, which needs no asymptotic distribution at all and is exact
under iid gaps. For m > 1 the statistic is the unweighted SUM of per-segment eq (7) values;
there is no eq (16) analogue to borrow, because eq (14)'s optimal weights are derived for a
normal limit these statistics do not have.

The bench PROMOTEd the permutation form (worst z = 2.67 against 3.42, mean power 0.612 at
n = 100) and REJECTed the asymptotic one.
