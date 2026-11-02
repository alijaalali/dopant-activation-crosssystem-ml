import numpy as np

EPS = 0.02


def logit_act(pct, eps=EPS):
    p = np.clip(np.asarray(pct, float) / 100.0, eps, 1 - eps)
    return np.log(p / (1 - p))


def inv_logit(z):
    return 1.0 / (1.0 + np.exp(-np.asarray(z, float)))


def gibbs_mixed(X, y, groups, n_iter=6000, burn=2000, seed=0, beta_sd=5.0, a0=2.0, b0=1.0,
                return_draws=False, disperse=False):
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
    if disperse:
        s2e = float(np.exp(rng.uniform(-2.0, 1.5)))
        s2u = [float(np.exp(rng.uniform(-3.0, 1.0))) for _ in K]
        u = [rng.normal(0.0, np.sqrt(v), size=k) for v, k in zip(s2u, K)]
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


def _autocov(x):
    n = len(x)
    x = np.asarray(x, float) - np.mean(x)
    f = np.fft.rfft(x, 2 * n)
    return np.fft.irfft(f * np.conj(f))[:n] / n * n / (n - 1)


def rhat_ess(chains):
    """Split R-hat and effective sample size for an array of shape (chains, draws)."""
    ch = np.asarray(chains, float)
    m, n = ch.shape
    half = n // 2
    sp = np.concatenate([ch[:, :half], ch[:, half:2 * half]], axis=0)
    nn = sp.shape[1]
    W = sp.var(axis=1, ddof=1).mean()
    B = nn * sp.mean(axis=1).var(ddof=1)
    var_plus = (nn - 1) / nn * W + B / nn
    rhat = float(np.sqrt(var_plus / W)) if W > 0 else float("nan")
    W2 = ch.var(axis=1, ddof=1).mean()
    B2 = n * ch.mean(axis=1).var(ddof=1) if m > 1 else 0.0
    vp = (n - 1) / n * W2 + B2 / n
    acov = np.mean([_autocov(c) for c in ch], axis=0)
    rho = 1.0 - (W2 - acov) / vp
    rho[0] = 1.0
    total, t = 0.0, 0
    while t + 1 < n:
        pair = rho[t] + rho[t + 1]
        if pair < 0:
            break
        total += pair
        t += 2
    tau = max(-1.0 + 2.0 * total, 1.0 / np.log10(max(m * n, 10)))
    return rhat, float(m * n / tau)


def load_sic_al(db_path):
    import pandas as pd
    db = pd.read_csv(db_path)
    al = db[(db.system == "SiC") & (db.dopant == "Al")].copy()
    al["log_time"] = np.log10(al.anneal_time_min)
    al["log_conc"] = np.log10(al.implanted_conc_cm3)
    need = ["activation_pct", "anneal_temp_C", "anneal_time_min", "implant_temp_C", "implanted_conc_cm3"]
    al = al.dropna(subset=need)
    al = al[(al.definition_class != "approximate") & (~al.two_step_anneal.fillna(False).astype(bool))]
    return al.reset_index(drop=True)


FEATS = ["implant_temp_C", "anneal_temp_C", "log_time", "log_conc"]


def _stat(v, ch=None):
    out = {"mean": float(np.mean(v)), "median": float(np.median(v)), "sd": float(np.std(v)),
           "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]}
    if ch is not None and ch.shape[0] > 1:
        out["rhat"], out["ess"] = rhat_ess(ch)
    return out


def fit_hier(frame, feats=FEATS, factors="both", beta_sd=5.0, eps=EPS, n_chains=4, n_iter=8000, burn=3000, seed0=1,
             a0=2.0, b0=1.0):
    """Hierarchical logit model fitted with several dispersed-start Gibbs chains.

    factors: 'both' (source and laboratory intercepts), 'source' or 'lab'.
    """
    X = frame[feats].to_numpy(float)
    mu, sd = X.mean(0), X.std(0, ddof=1)
    A = np.column_stack([np.ones(len(frame)), (X - mu) / sd])
    y = logit_act(frame.activation_pct, eps)
    src = pd_factorize(frame.source_id)
    lab = pd_factorize(frame.lab_group.fillna(frame.source_id))
    groups = {"both": [src, lab], "source": [src], "lab": [lab]}[factors]
    runs = [gibbs_mixed(A, y, groups, n_iter=n_iter, burn=burn, seed=seed0 + c, beta_sd=beta_sd, disperse=True, a0=a0, b0=b0)
            for c in range(n_chains)]
    nk = n_iter - burn
    beta = np.concatenate([r["beta"] for r in runs])
    s2e = np.concatenate([r["s2e"] for r in runs])
    s2u = np.concatenate([r["s2u"] for r in runs])
    names = ["intercept"] + list(feats)
    out = {"n_rows": int(len(frame)), "n_sources": int(src.max() + 1), "n_labs": int(lab.max() + 1),
           "n_chains": n_chains, "draws_per_chain": nk, "factors": factors, "eps": eps, "variance_prior_inv_gamma": [a0, b0]}
    out["coef"] = {}
    for j, nm in enumerate(names):
        ch = np.array([r["beta"][:, j] for r in runs])
        st = _stat(beta[:, j], ch)
        st["p_positive"] = float(np.mean(beta[:, j] > 0))
        out["coef"][nm] = st
    labels = {"both": ["source", "lab"], "source": ["source"], "lab": ["lab"]}[factors]
    out["sd"] = {}
    for k, nm in enumerate(labels):
        ch = np.array([np.sqrt(r["s2u"][:, k]) for r in runs])
        out["sd"][nm] = _stat(np.sqrt(s2u[:, k]), ch)
    out["sd"]["residual"] = _stat(np.sqrt(s2e), np.array([np.sqrt(r["s2e"]) for r in runs]))
    tot = s2u.sum(axis=1) + s2e
    out["share"] = {"between_groups_total": _stat(s2u.sum(axis=1) / tot), "residual": _stat(s2e / tot)}
    for k, nm in enumerate(labels):
        out["share"][nm] = _stat(s2u[:, k] / tot)
    if factors == "both":
        out["corr_source_lab_variance"] = float(np.corrcoef(s2u[:, 0], s2u[:, 1])[0, 1])
        out["lab_share_of_group_variance"] = _stat(s2u[:, 1] / s2u.sum(axis=1))
    pp = {}
    for j, f in enumerate(feats):
        a, b = np.zeros(A.shape[1]), np.zeros(A.shape[1])
        a[0] = b[0] = 1
        a[1 + j], b[1 + j] = 0.5, -0.5
        eff = 100 * (inv_logit(beta @ a) - inv_logit(beta @ b))
        pp[f] = {"mean": float(eff.mean()), "ci95": [float(np.percentile(eff, 2.5)), float(np.percentile(eff, 97.5))],
                 "p_positive": float(np.mean(eff > 0))}
    out["marginal_effect_pp_per_sd"] = pp
    out["feature_sd"] = {f: float(v) for f, v in zip(feats, sd)}
    out["_draws"] = {"beta": beta, "s2e": s2e, "s2u": s2u, "mu": mu, "sd": sd}
    return out


def pd_factorize(series):
    import pandas as pd
    return pd.factorize(series)[0]
