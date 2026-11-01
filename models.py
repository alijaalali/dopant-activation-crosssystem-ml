import numpy as np

EPS = 0.02


def logit_act(pct, eps=EPS):
    p = np.clip(np.asarray(pct, float) / 100.0, eps, 1 - eps)
    return np.log(p / (1 - p))


def inv_logit(z):
    return 1.0 / (1.0 + np.exp(-np.asarray(z, float)))


def gibbs_mixed(X, y, groups, n_iter=6000, burn=2000, seed=0, beta_sd=5.0, a0=2.0, b0=1.0,
                return_draws=False):
    """Gaussian linear mixed model with independent random intercepts for each grouping factor.

    X: (n, p) fixed-effect matrix including the intercept column.
    groups: list of integer-coded arrays (one per random-effect factor).
    Priors: beta ~ N(0, beta_sd^2); variances ~ InvGamma(a0, b0).
    """
    rng = np.random.default_rng(seed)
    n, p = X.shape
    G = [np.asarray(g) for g in groups]
    K = [int(g.max()) + 1 for g in G]
    beta = np.zeros(p)
    u = [np.zeros(k) for k in K]
    s2e = float(np.var(y)) + 1e-3
    s2u = [0.5 for _ in K]
    XtX = X.T @ X
    out_beta, out_s2e, out_s2u, out_u = [], [], [], []
    for it in range(n_iter):
        ru = np.zeros(n)
        for g, uu in zip(G, u):
            ru += uu[g]
        prec = XtX / s2e + np.eye(p) / beta_sd ** 2
        cov = np.linalg.inv(prec)
        mean = cov @ (X.T @ (y - ru)) / s2e
        beta = rng.multivariate_normal(mean, cov)
        fx = X @ beta
        for j, (g, k) in enumerate(zip(G, K)):
            other = np.zeros(n)
            for jj, (g2, uu2) in enumerate(zip(G, u)):
                if jj != j:
                    other += uu2[g2]
            r = y - fx - other
            cnt = np.bincount(g, minlength=k)
            sm = np.bincount(g, weights=r, minlength=k)
            v = 1.0 / (cnt / s2e + 1.0 / s2u[j])
            u[j] = rng.normal(v * sm / s2e, np.sqrt(v))
            s2u[j] = (b0 + 0.5 * np.sum(u[j] ** 2)) / rng.gamma(a0 + 0.5 * k)
        res = y - fx
        for g, uu in zip(G, u):
            res = res - uu[g]
        s2e = (b0 + 0.5 * np.sum(res ** 2)) / rng.gamma(a0 + 0.5 * n)
        if it >= burn:
            out_beta.append(beta.copy())
            out_s2e.append(s2e)
            out_s2u.append(list(s2u))
            if return_draws:
                out_u.append([uu.copy() for uu in u])
    res = {"beta": np.array(out_beta), "s2e": np.array(out_s2e), "s2u": np.array(out_s2u)}
    if return_draws:
        res["u"] = out_u
    return res


def summarize(draws, names):
    b = draws["beta"]
    return {nm: {"mean": float(b[:, j].mean()), "sd": float(b[:, j].std()),
                 "ci95": [float(np.percentile(b[:, j], 2.5)), float(np.percentile(b[:, j], 97.5))],
                 "p_positive": float(np.mean(b[:, j] > 0))} for j, nm in enumerate(names)}


def reml_profile(X, y, groups_int):
    """REML fit of a random-intercept model (one factor); returns beta, s2e, s2u and the REML log-likelihood."""
    from scipy.optimize import minimize_scalar
    n, p = X.shape
    g = np.asarray(groups_int)
    idx = [np.where(g == k)[0] for k in range(g.max() + 1)]

    def fit(logr):
        V = np.eye(n)
        for ii in idx:
            V[np.ix_(ii, ii)] += np.exp(logr)
        Vi = np.linalg.inv(V)
        M = X.T @ Vi @ X
        beta = np.linalg.solve(M, X.T @ Vi @ y)
        res = y - X @ beta
        s2 = float(res @ Vi @ res) / (n - p)
        ll = -0.5 * ((n - p) * np.log(s2) + np.linalg.slogdet(V)[1] + np.linalg.slogdet(M)[1])
        return beta, s2, np.exp(logr) * s2, ll, np.linalg.inv(M) * s2

    opt = minimize_scalar(lambda t: -fit(t)[3], bounds=(-8, 6), method="bounded")
    return fit(opt.x)


KB_EV = 8.617333262e-5


def arrhenius_log_time(T_C, t_min, Ea_eV, Tref_C=1700.0):
    T = np.asarray(T_C, float) + 273.15
    Tr = Tref_C + 273.15
    return np.log(np.asarray(t_min, float)) - (Ea_eV / KB_EV) * (1.0 / T - 1.0 / Tr)
