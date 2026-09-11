"""Generate the reference values used by tests/testthat/test-power_cc.R and
tests/testthat/test-power_nri.R.

Everything below is coded from the formulas printed in the article, not from the
R sources of this package.  The point of the two CSV files is to check the R
implementation against a different language, a different numerical library and a
different summation order, so this script deliberately does not import ssNRI or
translate its code.

Requires numpy and scipy.  Run from the package root:

    python inst/reference/generate_reference_values.py

It rewrites tests/testthat/ref_power_cc.csv and tests/testthat/ref_power_nri.csv
in place and prints the largest absolute change against the files already there.
"""

import csv
import os
import numpy as np
from scipy.stats import norm, binom

ALPHA_DEFAULT = 0.025
HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(HERE, os.pardir, os.pardir, "tests", "testthat")


def joint_cells(pi, omega, rho):
    """The four cells of Table 1: (no response, complete) and (response, complete)."""
    s = np.sqrt(pi * (1 - pi) * omega * (1 - omega))
    cell_00 = (1 - pi) * (1 - omega) + rho * s
    cell_10 = pi * (1 - omega) - rho * s
    return cell_00, cell_10


def power_cc(n1, n0, pi1, pi0, omega1, omega0, rho1, rho0, alpha):
    """Exact power and size of the one-sided pooled Z-test under complete case
    analysis, by enumeration over completer counts and responder counts.

    Under the null the response probability is replaced by the sample-size
    weighted pooled value; the dropout probabilities are design inputs and stay
    at their group-specific values.
    """
    c1_00, c1_10 = joint_cells(pi1, omega1, rho1)
    c0_00, c0_10 = joint_cells(pi0, omega0, rho0)
    pi1_cc = c1_10 / (c1_00 + c1_10)
    pi0_cc = c0_10 / (c0_00 + c0_10)
    pi_bar = (n1 * pi1_cc + n0 * pi0_cc) / (n1 + n0)

    z = norm.ppf(1 - alpha)
    completers1 = binom.pmf(np.arange(n1 + 1), n1, 1 - omega1)
    completers0 = binom.pmf(np.arange(n0 + 1), n0, 1 - omega0)

    power = 0.0
    size = 0.0
    for m1 in range(1, n1 + 1):
        x1 = np.arange(m1 + 1)
        f1_alt = binom.pmf(x1, m1, pi1_cc)
        f1_null = binom.pmf(x1, m1, pi_bar)
        for m0 in range(1, n0 + 1):
            x0 = np.arange(m0 + 1)
            grid1 = x1[:, None]
            grid0 = x0[None, :]
            pooled = (grid1 + grid0) / (m1 + m0)
            se = np.sqrt(pooled * (1 - pooled) * (1.0 / m1 + 1.0 / m0))
            se = np.where(se < 1e-10, 1e-10, se)
            reject = ((grid1 / m1 - grid0 / m0) / se) > z
            weight = completers1[m1] * completers0[m0]
            power += weight * float(
                (f1_alt[:, None] * binom.pmf(x0, m0, pi0_cc)[None, :] * reject).sum()
            )
            size += weight * float(
                (f1_null[:, None] * binom.pmf(x0, m0, pi_bar)[None, :] * reject).sum()
            )
    return pi1_cc, pi0_cc, power, size


def power_nri(n1, n0, pi1, pi0, omega1, omega0, rho1, rho0, alpha):
    """Exact power and size of the same test under non-responder imputation,
    where every randomized patient contributes and the denominators are fixed."""
    pi1_nri = joint_cells(pi1, omega1, rho1)[1]
    pi0_nri = joint_cells(pi0, omega0, rho0)[1]
    pi_bar = (n1 * pi1_nri + n0 * pi0_nri) / (n1 + n0)

    z = norm.ppf(1 - alpha)
    x1 = np.arange(n1 + 1)[:, None]
    x0 = np.arange(n0 + 1)[None, :]
    pooled = (x1 + x0) / (n1 + n0)
    se = np.sqrt(pooled * (1 - pooled) * (1.0 / n1 + 1.0 / n0))
    se = np.where(se < 1e-10, 1e-10, se)
    reject = ((x1 / n1 - x0 / n0) / se) > z

    f1_alt = binom.pmf(np.arange(n1 + 1), n1, pi1_nri)[:, None]
    f0_alt = binom.pmf(np.arange(n0 + 1), n0, pi0_nri)[None, :]
    f1_null = binom.pmf(np.arange(n1 + 1), n1, pi_bar)[:, None]
    f0_null = binom.pmf(np.arange(n0 + 1), n0, pi_bar)[None, :]
    power = float((f1_alt * f0_alt * reject).sum())
    size = float((f1_null * f0_null * reject).sum())
    return pi1_nri, pi0_nri, power, size


CC_CONFIGS = [
    # n1, n0, pi1, pi0, omega1, omega0, rho1, rho0, alpha
    (40, 40, 0.6, 0.4, 0.1, 0.1, 0, 0, 0.025),
    (50, 50, 0.35, 0.2, 0.05, 0.05, 0, 0, 0.025),
    (45, 45, 0.45, 0.3, 0.2, 0.2, 0.3, 0.3, 0.025),
    (45, 45, 0.45, 0.3, 0.2, 0.2, -0.3, -0.3, 0.025),
    (40, 60, 0.55, 0.35, 0.1, 0.25, 0.2, -0.1, 0.025),
    (30, 45, 0.5, 0.25, 0.15, 0.05, -0.1, 0.35, 0.05),
    (55, 55, 0.7, 0.5, 0.3, 0.1, 0, 0, 0.025),
    (25, 25, 0.4, 0.2, 0.02, 0.02, 0, 0, 0.025),
    (35, 35, 0.35, 0.2, 0.05, 0.05, 0.05, 0.05, 0.025),
    (60, 40, 0.65, 0.45, 0.12, 0.18, 0.15, -0.2, 0.01),
    (80, 80, 0.45, 0.25, 0.15, 0.15, 0, 0, 0.025),
    (70, 70, 0.5, 0.3, 0.1, 0.1, 0.2, 0.2, 0.025),
    (100, 90, 0.4, 0.25, 0.08, 0.2, -0.15, 0.1, 0.025),
]

NRI_CONFIGS = [
    # n1, n0, pi1, pi0, omega1, omega0, rho1, rho0, alpha
    (60, 60, 0.6, 0.4, 0.2, 0.2, 0, 0, 0.025),
    (60, 60, 0.4, 0.4, 0.05, 0.2, 0, 0, 0.025),
    (80, 80, 0.45, 0.3, 0.1, 0.1, 0.25, -0.2, 0.025),
    (120, 80, 0.35, 0.2, 0.15, 0.05, -0.15, 0.3, 0.025),
    (200, 200, 0.47, 0.41, 0.1, 0.1, 0, 0, 0.025),
    (55, 70, 0.55, 0.35, 0.25, 0.12, 0.1, 0.1, 0.01),
]


def write(path, configs, worker, label):
    head = ["n1", "n0", "pi1", "pi0", "omega1", "omega0", "rho1", "rho0", "alpha",
            "pi1_" + label, "pi0_" + label, "power", "type1_error"]
    rows = []
    for cfg in configs:
        n1, n0, pi1, pi0, om1, om0, rho1, rho0, alpha = cfg
        p1, p0, pw, t1 = worker(n1, n0, pi1, pi0, om1, om0, rho1, rho0, alpha)
        rows.append([n1, n0] + [repr(float(v)) for v in
                     (pi1, pi0, om1, om0, rho1, rho0, alpha, p1, p0, pw, t1)])
    old = {}
    if os.path.exists(path):
        with open(path, newline="") as handle:
            for i, rec in enumerate(csv.DictReader(handle)):
                old[i] = rec
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(head)
        writer.writerows(rows)
    worst = 0.0
    for i, row in enumerate(rows):
        if i in old:
            for key, new in zip(head[9:], row[9:]):
                worst = max(worst, abs(float(new) - float(old[i][key])))
    print("%s: %d rows, largest change against the previous file %.3e"
          % (os.path.basename(path), len(rows), worst))


if __name__ == "__main__":
    write(os.path.join(TESTS, "ref_power_cc.csv"), CC_CONFIGS, power_cc, "CC")
    write(os.path.join(TESTS, "ref_power_nri.csv"), NRI_CONFIGS, power_nri, "NRI")
