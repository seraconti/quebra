"""Rigdon and Basu's failure-truncated example, as a test fixture.

Source: S. E. Rigdon and A. P. Basu, Statistical Methods for the Reliability of Repairable
Systems, Wiley (2000).

- `FAILURE_TIMES`: Table 4.1, p. 88, system 3. Twelve cumulative failure times, as printed,
  failure truncated at the twelfth, `t_n = 227`.
- `PUBLISHED`: Example 4.8, pp. 110-111. The Laplace statistic for failure-truncated data,
  `L = 2.23`, and its `P = 0.0257`. The example's own hand check: the first 11 times sum to
  1734, and `(1734/11 - 227/2) / (227 / sqrt(12 * 11)) = 2.23`.

Failure truncation is the case `checks.result.Segment.ends_at_event` names: observation
stops at the last event, so the other `n - 1` event times are what the statistic reads.
"""

from __future__ import annotations

import numpy as np

FAILURE_TIMES = np.array(
    [45, 76, 113, 129, 152, 174, 193, 199, 210, 219, 224, 227], dtype=float
)

PUBLISHED = {
    "n": 12,
    "t_n": 227.0,
    "sum_first_n_minus_1": 1734.0,
    "laplace": 2.23,
    "p_value": 0.0257,
}
