#!/usr/bin/env python3
"""Fisher and profile-likelihood sketches for a shared metabolic ODE.

Two observation schedules are compared on the same right-hand side:

  S_ratio     steady lactate/glucose ratio (single summary)
  M_snapshot  four steady channels, correlated noise, lineage loading

A third schedule, D_transient, is an upper bound only. It is a short
time course with independent noise. It is not a CCLE-style snapshot.

The metabolite draws are a synthetic surrogate. They are not a download
of CCLE or DepMap files.

Research computation only. Not a medical device, dose, or clinical tool.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)

SEED = 20260921
THETA_NAMES = ["k_gp", "k_pl", "k_pq", "d_l", "d_q"]
THETA0 = np.array([0.80, 0.50, 0.30, 0.40, 0.60], dtype=float)
CHANNELS = ["G", "P", "L", "Q"]
LINEAGES = ["BRCA", "LUAD", "COAD", "AML"]
U_G = {"BRCA": 1.00, "LUAD": 1.40, "COAD": 0.80, "AML": 1.20}
U_Q = {"BRCA": 0.40, "LUAD": 0.55, "COAD": 0.50, "AML": 0.30}
N_REP = 6
N_KIN = THETA0.size

# Synthetic correlation on the log-concentration scale. Not estimated from CCLE.
CORR = np.array(
    [
        [1.00, 0.55, 0.45, 0.15],
        [0.55, 1.00, 0.60, 0.20],
        [0.45, 0.60, 1.00, 0.10],
        [0.15, 0.20, 0.10, 1.00],
    ],
    dtype=float,
)
SIGMA_LOG = 0.15
SIGMA_RATIO = 0.08
SIGMA_BATCH = 0.20
Y0 = np.array([1.20, 0.40, 0.30, 0.50], dtype=float)
T_END = 15.0
CHI2_95 = 3.841
RANK_TOL = 1e-8
PRACTICAL_EIG_FRAC = 1e-3

BOUNDS_LO = np.array([0.05, 0.05, 0.02, 0.05, 0.05], dtype=float)
BOUNDS_HI = np.array([4.0, 4.0, 3.0, 4.0, 4.0], dtype=float)
BETA_LO = -1.5
BETA_HI = 1.5


def steady(theta: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    k_gp, k_pl, k_pq, d_l, d_q = (float(v) for v in theta)
    if min(k_gp, k_pl, k_pq, d_l, d_q) <= 0.0:
        return np.full(4, np.nan)
    g = u_g / k_gp
    p = u_g / (k_pl + k_pq)
    lac = k_pl * p / d_l
    q = (u_q + k_pq * p) / d_q
    out = np.array([g, p, lac, q], dtype=float)
    if np.any(out <= 0.0) or not np.all(np.isfinite(out)):
        return np.full(4, np.nan)
    return out


def ratio(theta: np.ndarray, u_g: float, u_q: float) -> float:
    x = steady(theta, u_g, u_q)
    if not np.all(np.isfinite(x)):
        return float("nan")
    return float(x[2] / x[0])


def rhs(t: float, y: np.ndarray, theta: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    g, p, lac, q = y
    k_gp, k_pl, k_pq, d_l, d_q = theta
    return np.array(
        [
            u_g - k_gp * g,
            k_gp * g - (k_pl + k_pq) * p,
            k_pl * p - d_l * lac,
            u_q + k_pq * p - d_q * q,
        ],
        dtype=float,
    )


def simulate(theta: np.ndarray, u_g: float, u_q: float, t_eval: np.ndarray) -> np.ndarray:
    sol = solve_ivp(
        lambda t, y: rhs(t, y, theta, u_g, u_q),
        (float(t_eval[0]), float(t_eval[-1])),
        Y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
        max_step=0.25,
    )
    if (not sol.success) or sol.y.shape[1] != t_eval.size or not np.all(np.isfinite(sol.y)):
        return np.full((t_eval.size, 4), np.nan)
    return sol.y.T


def sigma_matrix() -> np.ndarray:
    return (SIGMA_LOG**2) * CORR


def whiten_matrix(cov: np.ndarray) -> np.ndarray:
    chol = cholesky(cov, lower=True)
    eye = np.eye(cov.shape[0])
    return solve_triangular(chol, eye, lower=True)


def eig_report(mat: np.ndarray) -> dict:
    evals = np.sort(np.linalg.eigvalsh(0.5 * (mat + mat.T)))[::-1]
    evals = np.clip(evals, 0.0, None)
    smax = float(evals[0]) if evals.size and evals[0] > 0 else 0.0
    tol = RANK_TOL * smax if smax > 0 else 1e-18
    rank = int(np.sum(evals > tol))
    practical_tol = PRACTICAL_EIG_FRAC * smax if smax > 0 else 1e-18
    practical_rank = int(np.sum(evals > practical_tol))
    if rank >= 1 and evals[rank - 1] > 0:
        cond = float(evals[0] / evals[rank - 1])
    else:
        cond = float("inf")
    return {
        "rank": rank,
        "practical_rank": practical_rank,
        "n_params": int(mat.shape[0]),
        "eigenvalues": [float(v) for v in evals],
        "log10_eigenvalues": [float(np.log10(v + 1e-30)) for v in evals],
        "condition_identifiable_block": cond,
    }


def kinetic_jacobian_log(theta: np.ndarray, u_g: float, u_q: float) -> np.ndarray:
    """Finite-difference ∂ log x* / ∂θ, shape (4, n_kin)."""
    x0 = steady(theta, u_g, u_q)
    jac = np.zeros((4, N_KIN), dtype=float)
    if not np.all(np.isfinite(x0)):
        return jac
    for j in range(N_KIN):
        step = 1e-6 * max(abs(float(theta[j])), 1e-3)
        th = theta.copy()
        th[j] += step
        xp = steady(th, u_g, u_q)
        if not np.all(np.isfinite(xp)):
            jac[:, j] = np.nan
            continue
        jac[:, j] = (np.log(xp) - np.log(x0)) / step
    return jac


def fim_ratio(theta: np.ndarray) -> dict:
    """Scalar log-ratio, independent Gaussian noise, no lineage loading."""
    f = np.zeros((N_KIN, N_KIN), dtype=float)
    grads = {}
    for name in LINEAGES:
        s0 = ratio(theta, U_G[name], U_Q[name])
        g = np.zeros(N_KIN, dtype=float)
        for j in range(N_KIN):
            step = 1e-6 * max(abs(float(theta[j])), 1e-3)
            th = theta.copy()
            th[j] += step
            sp = ratio(th, U_G[name], U_Q[name])
            g[j] = (np.log(sp) - np.log(s0)) / step
        grads[name] = [float(v) for v in g]
        f += N_REP * np.outer(g, g) / (SIGMA_RATIO**2)
    rep = eig_report(f)
    rep.update(
        {
            "schedule": "S_ratio",
            "parameter_block": "kinetics",
            "n_lineages": len(LINEAGES),
            "n_replicates": N_REP,
            "n_observations": len(LINEAGES) * N_REP,
            "gradient_log_ratio": grads,
            "zero_columns": [
                THETA_NAMES[j]
                for j in range(N_KIN)
                if float(np.linalg.norm(f[:, j])) <= 1e-8 * (np.linalg.norm(f) + 1.0)
            ],
        }
    )
    if rep["rank"] == N_KIN:
        cov = np.linalg.inv(f)
        rep["crlb_sd"] = np.sqrt(np.clip(np.diag(cov), 0, None)).tolist()
        rep["crlb_rel"] = (np.array(rep["crlb_sd"]) / THETA0).tolist()
    else:
        rep["crlb_sd"] = None
        rep["crlb_rel"] = None
    return rep


def fim_snapshot(theta: np.ndarray) -> dict:
    """Four log-channels. Lineage loading is a fixed nuisance, then profiled out."""
    n_beta = len(LINEAGES)
    n_p = N_KIN + n_beta
    f = np.zeros((n_p, n_p), dtype=float)
    w = whiten_matrix(sigma_matrix())
    for i, name in enumerate(LINEAGES):
        jk = kinetic_jacobian_log(theta, U_G[name], U_Q[name])
        j = np.zeros((4, n_p), dtype=float)
        j[:, :N_KIN] = jk
        j[:, N_KIN + i] = 1.0
        jw = w @ j
        f += N_REP * (jw.T @ jw)
    full = eig_report(f)
    ftt = f[:N_KIN, :N_KIN]
    ftb = f[:N_KIN, N_KIN:]
    fbb = f[N_KIN:, N_KIN:]
    schur = ftt - ftb @ np.linalg.solve(fbb, ftb.T)
    prof = eig_report(schur)
    # Smallest profiled eigenvector: which kinetics are weakest.
    evals, evecs = np.linalg.eigh(0.5 * (schur + schur.T))
    order = np.argsort(evals)
    weak = evecs[:, order[0]]
    weak = weak / (np.linalg.norm(weak) + 1e-15)
    out = {
        "schedule": "M_snapshot",
        "full_with_batch": full,
        "profiled_kinetics": prof,
        "weakest_profiled_direction": {
            THETA_NAMES[j]: float(weak[j]) for j in range(N_KIN)
        },
        "n_lineages": len(LINEAGES),
        "n_replicates": N_REP,
        "n_channels": 4,
        "n_batch": n_beta,
        "n_observations": len(LINEAGES) * N_REP * 4,
    }
    if prof["rank"] == N_KIN:
        cov = np.linalg.inv(schur)
        out["crlb_sd"] = np.sqrt(np.clip(np.diag(cov), 0, None)).tolist()
        out["crlb_rel"] = (np.array(out["crlb_sd"]) / THETA0).tolist()
    else:
        out["crlb_sd"] = None
        out["crlb_rel"] = None
    return out


def fim_transient(theta: np.ndarray) -> dict:
    """Upper bound: BRCA time course, all four states, independent noise, no batch."""
    t_eval = np.linspace(0.0, T_END, 31)
    name = "BRCA"
    y0 = simulate(theta, U_G[name], U_Q[name], t_eval)
    sens = np.zeros((t_eval.size, 4, N_KIN), dtype=float)
    for j in range(N_KIN):
        step = 1e-5 * max(abs(float(theta[j])), 1e-3)
        th = theta.copy()
        th[j] += step
        yp = simulate(th, U_G[name], U_Q[name], t_eval)
        sens[:, :, j] = (yp - y0) / step
    f = np.zeros((N_KIN, N_KIN), dtype=float)
    for k in range(4):
        sig = 0.05 * np.abs(y0[:, k]) + 0.02
        w = 1.0 / sig**2
        sj = sens[:, k, :]
        f += (sj * w[:, None]).T @ sj
    rep = eig_report(f)
    rep.update(
        {
            "schedule": "D_transient",
            "lineage": name,
            "n_times": int(t_eval.size),
            "channels": CHANNELS,
            "noise": "independent 5% CV plus 0.02 floor",
            "role": "upper bound, not a CCLE-style snapshot",
        }
    )
    if rep["rank"] == N_KIN:
        cov = np.linalg.inv(f)
        rep["crlb_sd"] = np.sqrt(np.clip(np.diag(cov), 0, None)).tolist()
        rep["crlb_rel"] = (np.array(rep["crlb_sd"]) / THETA0).tolist()
    else:
        rep["crlb_sd"] = None
        rep["crlb_rel"] = None
    return rep


def _pack(free: np.ndarray, fixed_index: int, fixed_value: float, with_beta: bool) -> tuple[np.ndarray, np.ndarray | None]:
    n_beta = len(LINEAGES) if with_beta else 0
    theta = np.empty(N_KIN, dtype=float)
    mask = np.ones(N_KIN, dtype=bool)
    mask[fixed_index] = False
    theta[mask] = free[: N_KIN - 1]
    theta[fixed_index] = fixed_value
    beta = free[N_KIN - 1 :] if with_beta else None
    if with_beta and beta.size != n_beta:
        raise RuntimeError("beta length")
    return theta, beta


def residuals_ratio(theta: np.ndarray, y_log: np.ndarray) -> np.ndarray:
    chunks = []
    for i, name in enumerate(LINEAGES):
        s = ratio(theta, U_G[name], U_Q[name])
        if not np.isfinite(s) or s <= 0:
            return np.full(y_log.size, 1e3)
        pred = np.log(s)
        chunks.append((y_log[i] - pred) / SIGMA_RATIO)
    r = np.concatenate(chunks)
    if not np.all(np.isfinite(r)):
        return np.full(y_log.size, 1e3)
    return r


def residuals_snapshot(theta: np.ndarray, beta: np.ndarray, y_log: np.ndarray, w: np.ndarray) -> np.ndarray:
    chunks = []
    for i, name in enumerate(LINEAGES):
        x = steady(theta, U_G[name], U_Q[name])
        if not np.all(np.isfinite(x)):
            return np.full(y_log.size, 1e3)
        pred = np.log(x) + beta[i]
        for r in range(y_log.shape[1]):
            chunks.append(w @ (y_log[i, r] - pred))
    out = np.concatenate(chunks)
    if not np.all(np.isfinite(out)):
        return np.full(out.size, 1e3)
    return out


def _call_profile(grid: np.ndarray, chi2: np.ndarray) -> dict:
    chi2 = np.asarray(chi2, dtype=float)
    chi2_min = float(np.nanmin(chi2))
    thresh = chi2_min + CHI2_95
    finite = chi2 <= thresh
    if np.any(finite):
        lo = float(grid[finite].min())
        hi = float(grid[finite].max())
    else:
        lo, hi = float("nan"), float("nan")
    bounded = bool(finite.size) and (not bool(finite[0])) and (not bool(finite[-1]))
    spread = float(np.nanmax(chi2) - np.nanmin(chi2))
    if spread < 0.5:
        call = "structurally_flat"
        bounded = False
    elif bounded:
        call = "identifiable_on_grid"
    else:
        call = "practical_nonidentifiable"
    return {
        "chi2_min": chi2_min,
        "threshold_95": thresh,
        "delta_chi2": [float(v - chi2_min) for v in chi2],
        "interval_95_on_grid": [lo, hi],
        "bounded_95": bounded,
        "chi2_spread": spread,
        "call": call,
    }


def profile_ratio(param_index: int, y_log: np.ndarray, n_grid: int = 13) -> dict:
    true = float(THETA0[param_index])
    grid = true * np.logspace(np.log10(0.45), np.log10(2.2), n_grid)
    chi2 = []
    for val in grid:
        def fun(free, val=val):
            theta, _ = _pack(free, param_index, float(val), with_beta=False)
            return residuals_ratio(theta, y_log)

        free0 = np.delete(THETA0, param_index)
        lo = np.delete(BOUNDS_LO, param_index)
        hi = np.delete(BOUNDS_HI, param_index)
        free0 = np.clip(free0, lo, hi)
        try:
            opt = least_squares(fun, free0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=250)
            chi2.append(float(opt.fun @ opt.fun))
        except ValueError:
            chi2.append(1e12)
    chi2_arr = np.asarray(chi2, dtype=float)
    out = {
        "schedule": "S_ratio",
        "parameter": THETA_NAMES[param_index],
        "true": true,
        "grid": [float(v) for v in grid],
        "chi2": [float(v) for v in chi2_arr],
    }
    out.update(_call_profile(grid, chi2_arr))
    return out


def fim_snapshot_gauged(theta: np.ndarray, gauge_index: int = 0) -> dict:
    """Kinetic FIM with one rate held fixed, so the scale symmetry is removed."""
    free = [j for j in range(N_KIN) if j != gauge_index]
    n_beta = len(LINEAGES)
    n_p = len(free) + n_beta
    f = np.zeros((n_p, n_p), dtype=float)
    w = whiten_matrix(sigma_matrix())
    for i, name in enumerate(LINEAGES):
        jk = kinetic_jacobian_log(theta, U_G[name], U_Q[name])
        j = np.zeros((4, n_p), dtype=float)
        for col, jpar in enumerate(free):
            j[:, col] = jk[:, jpar]
        j[:, len(free) + i] = 1.0
        jw = w @ j
        f += N_REP * (jw.T @ jw)
    ftt = f[: len(free), : len(free)]
    ftb = f[: len(free), len(free) :]
    fbb = f[len(free) :, len(free) :]
    schur = ftt - ftb @ np.linalg.solve(fbb, ftb.T)
    prof = eig_report(schur)
    names = [THETA_NAMES[j] for j in free]
    out = {
        "schedule": "M_snapshot_gauge",
        "gauge": THETA_NAMES[gauge_index],
        "gauge_value": float(theta[gauge_index]),
        "free_names": names,
        "profiled_kinetics": prof,
    }
    if prof["rank"] == len(free):
        cov = np.linalg.inv(schur)
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))
        out["crlb_sd"] = sd.tolist()
        out["crlb_rel"] = (sd / theta[free]).tolist()
    else:
        out["crlb_sd"] = None
        out["crlb_rel"] = None
    return out


def profile_snapshot_gauged(
    param_index: int,
    y_log: np.ndarray,
    w: np.ndarray,
    gauge_index: int = 0,
    n_grid: int = 21,
) -> dict:
    """Profile one kinetic parameter with the scale gauge held at THETA0."""
    if param_index == gauge_index:
        raise ValueError("gauge and profiled index coincide")
    true = float(THETA0[param_index])
    grid = true * np.logspace(np.log10(0.15), np.log10(3.0), n_grid)
    n_beta = len(LINEAGES)
    free_idx = [j for j in range(N_KIN) if j not in (param_index, gauge_index)]
    chi2 = []
    for val in grid:
        def fun(free, val=val):
            theta = THETA0.copy()
            theta[gauge_index] = THETA0[gauge_index]
            theta[param_index] = float(val)
            theta[free_idx] = free[: len(free_idx)]
            beta = free[len(free_idx) :]
            return residuals_snapshot(theta, beta, y_log, w)

        free0 = np.concatenate([THETA0[free_idx], np.zeros(n_beta)])
        lo = np.concatenate([BOUNDS_LO[free_idx], np.full(n_beta, BETA_LO)])
        hi = np.concatenate([BOUNDS_HI[free_idx], np.full(n_beta, BETA_HI)])
        free0 = np.clip(free0, lo, hi)
        try:
            opt = least_squares(fun, free0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=400)
            chi2.append(float(opt.fun @ opt.fun))
        except ValueError:
            chi2.append(1e12)
    chi2_arr = np.asarray(chi2, dtype=float)
    out = {
        "schedule": "M_snapshot_gauge",
        "gauge": THETA_NAMES[gauge_index],
        "parameter": THETA_NAMES[param_index],
        "true": true,
        "grid": [float(v) for v in grid],
        "chi2": [float(v) for v in chi2_arr],
    }
    out.update(_call_profile(grid, chi2_arr))
    return out


def profile_snapshot(param_index: int, y_log: np.ndarray, w: np.ndarray, n_grid: int = 13) -> dict:
    true = float(THETA0[param_index])
    grid = true * np.logspace(np.log10(0.45), np.log10(2.2), n_grid)
    n_beta = len(LINEAGES)
    chi2 = []
    for val in grid:
        def fun(free, val=val):
            theta, beta = _pack(free, param_index, float(val), with_beta=True)
            return residuals_snapshot(theta, beta, y_log, w)

        free0 = np.concatenate([np.delete(THETA0, param_index), np.zeros(n_beta)])
        lo = np.concatenate([np.delete(BOUNDS_LO, param_index), np.full(n_beta, BETA_LO)])
        hi = np.concatenate([np.delete(BOUNDS_HI, param_index), np.full(n_beta, BETA_HI)])
        free0 = np.clip(free0, lo, hi)
        try:
            opt = least_squares(fun, free0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=400)
            chi2.append(float(opt.fun @ opt.fun))
        except ValueError:
            chi2.append(1e12)
    chi2_arr = np.asarray(chi2, dtype=float)
    out = {
        "schedule": "M_snapshot",
        "parameter": THETA_NAMES[param_index],
        "true": true,
        "grid": [float(v) for v in grid],
        "chi2": [float(v) for v in chi2_arr],
    }
    out.update(_call_profile(grid, chi2_arr))
    return out


def slice_fixed(param_index: int, schedule: str, y_builder) -> dict:
    """One-at-a-time slice with every other parameter held at THETA0. Not a profile."""
    true = float(THETA0[param_index])
    grid = true * np.logspace(np.log10(0.45), np.log10(2.2), 41)
    chi2 = []
    for val in grid:
        theta = THETA0.copy()
        theta[param_index] = val
        r = y_builder(theta)
        chi2.append(float(r @ r))
    chi2_arr = np.asarray(chi2, dtype=float)
    out = {
        "schedule": schedule,
        "parameter": THETA_NAMES[param_index],
        "true": true,
        "kind": "slice_others_fixed",
        "grid": [float(v) for v in grid],
        "chi2": [float(v) for v in chi2_arr],
    }
    out.update(_call_profile(grid, chi2_arr))
    return out


def make_ratio_data(rng: np.random.Generator, noise: bool) -> np.ndarray:
    y = np.zeros((len(LINEAGES), N_REP), dtype=float)
    for i, name in enumerate(LINEAGES):
        mu = np.log(ratio(THETA0, U_G[name], U_Q[name]))
        if noise:
            y[i] = mu + SIGMA_RATIO * rng.normal(size=N_REP)
        else:
            y[i] = mu
    return y


def make_snapshot_data(rng: np.random.Generator, noise: bool) -> tuple[np.ndarray, np.ndarray]:
    cov = sigma_matrix()
    y = np.zeros((len(LINEAGES), N_REP, 4), dtype=float)
    beta = np.zeros(len(LINEAGES), dtype=float)
    for i, name in enumerate(LINEAGES):
        mu = np.log(steady(THETA0, U_G[name], U_Q[name]))
        if noise:
            beta[i] = SIGMA_BATCH * rng.normal()
            draws = rng.multivariate_normal(np.zeros(4), cov, size=N_REP)
            y[i] = mu + beta[i] + draws
        else:
            y[i] = mu
    return y, beta


def steady_table(theta: np.ndarray) -> dict:
    rows = {}
    for name in LINEAGES:
        x = steady(theta, U_G[name], U_Q[name])
        rows[name] = {
            "u_G": U_G[name],
            "u_Q": U_Q[name],
            "G": float(x[0]),
            "P": float(x[1]),
            "L": float(x[2]),
            "Q": float(x[3]),
            "L_over_G": float(x[2] / x[0]),
        }
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    corr_eigs = np.linalg.eigvalsh(CORR)
    if np.min(corr_eigs) <= 0:
        raise RuntimeError(f"CORR is not positive definite: {corr_eigs}")

    rows = steady_table(THETA0)
    ratios = [rows[n]["L_over_G"] for n in LINEAGES]
    fim_s = fim_ratio(THETA0)
    fim_m = fim_snapshot(THETA0)
    fim_d = fim_transient(THETA0)

    y_s0 = make_ratio_data(rng, noise=False)
    y_s1 = make_ratio_data(rng, noise=True)
    y_m0, beta0 = make_snapshot_data(rng, noise=False)
    y_m1, beta1 = make_snapshot_data(rng, noise=True)
    w = whiten_matrix(sigma_matrix())

    profiles = {}
    for pname in ("k_gp", "k_pl", "d_q"):
        idx = THETA_NAMES.index(pname)
        profiles[f"S|{pname}|noiseless"] = profile_ratio(idx, y_s0)
        profiles[f"S|{pname}|noisy"] = profile_ratio(idx, y_s1)
        profiles[f"M|{pname}|noiseless"] = profile_snapshot(idx, y_m0, w)
        profiles[f"M|{pname}|noisy"] = profile_snapshot(idx, y_m1, w)
    for pname in ("k_pl", "d_q"):
        idx = THETA_NAMES.index(pname)
        profiles[f"Mgauge|{pname}|noiseless"] = profile_snapshot_gauged(idx, y_m0, w)
        profiles[f"Mgauge|{pname}|noisy"] = profile_snapshot_gauged(idx, y_m1, w)
    fim_g = fim_snapshot_gauged(THETA0, gauge_index=0)

    def ratio_resid(theta: np.ndarray) -> np.ndarray:
        return residuals_ratio(theta, y_s0)

    def snap_resid(theta: np.ndarray) -> np.ndarray:
        return residuals_snapshot(theta, np.zeros(len(LINEAGES)), y_m0, w)

    slices = {
        "S|k_gp|slice": slice_fixed(0, "S_ratio", ratio_resid),
        "M|k_gp|slice": slice_fixed(0, "M_snapshot", snap_resid),
        "S|d_q|slice": slice_fixed(4, "S_ratio", ratio_resid),
        "M|d_q|slice": slice_fixed(4, "M_snapshot", snap_resid),
    }

    t_plot = np.linspace(0.0, T_END, 301)
    trajectories = {}
    for name in ("BRCA", "LUAD"):
        y = simulate(THETA0, U_G[name], U_Q[name], t_plot)
        trajectories[name] = {
            "t": t_plot[::5].tolist(),
            "y": y[::5].tolist(),
            "terminal": y[-1].tolist(),
            "steady": [rows[name][c] for c in CHANNELS],
        }

    out = {
        "seed": SEED,
        "data": "synthetic multi-channel surrogate, not a CCLE download",
        "theta_names": THETA_NAMES,
        "theta0": THETA0.tolist(),
        "lineages": LINEAGES,
        "n_replicates": N_REP,
        "sigma_log": SIGMA_LOG,
        "sigma_ratio": SIGMA_RATIO,
        "sigma_batch": SIGMA_BATCH,
        "correlation": CORR.tolist(),
        "correlation_eigenvalues": [float(v) for v in corr_eigs],
        "steady": rows,
        "ratio_identical_across_lineages": bool(np.allclose(ratios, ratios[0])),
        "ratio_value": float(ratios[0]),
        "beta_noisy": {LINEAGES[i]: float(beta1[i]) for i in range(len(LINEAGES))},
        "beta_noiseless": {LINEAGES[i]: float(beta0[i]) for i in range(len(LINEAGES))},
        "fim_ratio": fim_s,
        "fim_snapshot": fim_m,
        "fim_snapshot_gauge_kgp": fim_g,
        "fim_transient": fim_d,
        "profiles": profiles,
        "slices": slices,
        "disclaimer": "Computational research. Not a medical device, CDS, dose, or cure.",
    }
    # Trajectories are large; keep a compact copy for the figure routine only.
    (ROOT / "results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    _make_figures(t_plot, trajectories, rows, fim_s, fim_m, fim_d, profiles, slices)
    _print_report(out)
    print(f"wrote {ROOT / 'results.json'}")


def _print_report(out: dict) -> None:
    print("ratio", out["ratio_value"], "identical", out["ratio_identical_across_lineages"])
    print("steady")
    for name, row in out["steady"].items():
        print(f"  {name}: G={row['G']:.4f} P={row['P']:.4f} L={row['L']:.4f} Q={row['Q']:.4f}")
    fr = out["fim_ratio"]
    print(
        f"S rank {fr['rank']}/{fr['n_params']} practical {fr['practical_rank']} "
        f"cond {fr['condition_identifiable_block']:.6g} zero {fr['zero_columns']}"
    )
    print("S evals", ["%.4g" % v for v in fr["eigenvalues"]])
    fm = out["fim_snapshot"]
    pk = fm["profiled_kinetics"]
    print(
        f"M full rank {fm['full_with_batch']['rank']}/{fm['full_with_batch']['n_params']} "
        f"profiled {pk['rank']}/{pk['n_params']} practical {pk['practical_rank']} "
        f"cond {pk['condition_identifiable_block']:.6g}"
    )
    print("M profiled evals", ["%.4g" % v for v in pk["eigenvalues"]])
    print("M crlb_rel", None if fm["crlb_rel"] is None else ["%.4g" % v for v in fm["crlb_rel"]])
    print("weak", fm["weakest_profiled_direction"])
    fd = out["fim_transient"]
    print(
        f"D rank {fd['rank']}/{fd['n_params']} practical {fd['practical_rank']} "
        f"cond {fd['condition_identifiable_block']:.6g}"
    )
    print("D crlb_rel", None if fd["crlb_rel"] is None else ["%.4g" % v for v in fd["crlb_rel"]])
    fg = out["fim_snapshot_gauge_kgp"]
    pk = fg["profiled_kinetics"]
    print(
        f"M gauge {fg['gauge']}={fg['gauge_value']} rank {pk['rank']}/{pk['n_params']} "
        f"practical {pk['practical_rank']} cond {pk['condition_identifiable_block']:.6g}"
    )
    print("M gauge names", fg["free_names"])
    print("M gauge crlb_rel", None if fg["crlb_rel"] is None else ["%.4g" % v for v in fg["crlb_rel"]])
    print("M gauge evals", ["%.4g" % v for v in pk["eigenvalues"]])
    for key, pr in out["profiles"].items():
        print(
            f"profile {key}: {pr['call']} spread {pr['chi2_spread']:.4g} "
            f"bounded {pr['bounded_95']} interval {pr['interval_95_on_grid']}"
        )
    for key, pr in out["slices"].items():
        print(f"slice {key}: {pr['call']} spread {pr['chi2_spread']:.4g} bounded {pr['bounded_95']}")


def _make_figures(t_plot, trajectories, rows, fim_s, fim_m, fim_d, profiles, slices) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), sharex=True)
    colors = {"BRCA": "#1f4e79", "LUAD": "#b85c38"}
    for ax, label in zip(axes.ravel(), ["glucose G", "pyruvate P", "lactate L", "glutamine pool Q"]):
        ax.set_ylabel(label)
    # Resimulate at full resolution for the figure (same seed-free ODE).
    for name, color in colors.items():
        y = simulate(THETA0, U_G[name], U_Q[name], t_plot)
        for ax, j in zip(axes.ravel(), range(4)):
            ax.plot(t_plot, y[:, j], color=color, lw=1.6, label=name)
            ax.axhline(rows[name][CHANNELS[j]], color=color, lw=0.7, ls=":")
    for ax in axes.ravel():
        ax.legend(fontsize=8, frameon=False)
    axes[1, 0].set_xlabel("time")
    axes[1, 1].set_xlabel("time")
    fig.suptitle("Relaxation to the lineage steady state (synthetic inputs)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGDIR / "trajectories_steady.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    x = np.arange(len(LINEAGES))
    width = 0.18
    for j, (ch, color) in enumerate(zip(CHANNELS, ["#1f4e79", "#4c78a8", "#b85c38", "#5b8c5a"])):
        vals = [rows[n][ch] for n in LINEAGES]
        ax.bar(x + (j - 1.5) * width, vals, width=width, label=ch, color=color)
    ax.set_xticks(x, LINEAGES)
    ax.set_ylabel("steady concentration (model units)")
    ax.set_title("Noise-free steady state by synthetic lineage")
    ax.legend(ncol=4, frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "steady_by_lineage.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    series = [
        (np.array(fim_s["eigenvalues"]), "S ratio, kinetics"),
        (np.array(fim_m["profiled_kinetics"]["eigenvalues"]), "M snapshot, kinetics after batch"),
        (np.array(fim_d["eigenvalues"]), "D transient upper bound"),
    ]
    for ev, lab in series:
        ax.semilogy(np.arange(1, ev.size + 1), ev + 1e-30, marker="o", lw=1.5, label=lab)
    ax.set_xlabel("eigenvalue index (largest first)")
    ax.set_ylabel("Fisher eigenvalue")
    ax.set_title("Fisher spectra on the shared kinetic parameters")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fim_eigenspectra.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.8), sharey=True)
    panels = [
        (axes[0], [("S|k_gp|noiseless", "S profile"), ("M|k_gp|noiseless", "M profile"), ("S|k_gp|slice", "S slice, others fixed")], "k_gp"),
        (axes[1], [("S|d_q|noiseless", "S profile"), ("M|d_q|noiseless", "M profile"), ("S|d_q|slice", "S slice, others fixed")], "d_q"),
    ]
    store = {**profiles, **slices}
    for ax, keys, title in panels:
        for key, lab in keys:
            pr = store[key]
            ax.plot(pr["grid"], pr["delta_chi2"], lw=1.6, label=lab)
        ax.axhline(CHI2_95, color="0.4", ls=":", lw=1.0)
        ax.set_xlabel(title)
        ax.set_title(title)
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("delta chi-square")
    fig.suptitle("Profiles refit the other parameters. Slices do not.", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGDIR / "profiles_kgp_dq.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for key, lab in [
        ("M|k_pl|noiseless", "M k_pl, scale free"),
        ("M|d_q|noiseless", "M d_q, scale free"),
        ("Mgauge|k_pl|noiseless", "M k_pl, k_gp fixed"),
        ("Mgauge|d_q|noiseless", "M d_q, k_gp fixed"),
        ("Mgauge|k_pl|noisy", "M k_pl, k_gp fixed, noisy"),
        ("Mgauge|d_q|noisy", "M d_q, k_gp fixed, noisy"),
    ]:
        pr = profiles[key]
        ax.plot(np.array(pr["grid"]) / pr["true"], pr["delta_chi2"], lw=1.5, label=lab)
    ax.axhline(CHI2_95, color="0.4", ls=":", lw=1.0)
    ax.set_xlabel("parameter / true value")
    ax.set_ylabel("delta chi-square")
    ax.set_title("Scale-free profiles are flat. A k_gp gauge restores a peak.")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR / "profiles_snapshot_noise.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
