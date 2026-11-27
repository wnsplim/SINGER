#!/usr/bin/env python3

import argparse
import sys

import numpy as np


def scalar_ess(x):
    x = np.asarray(x, dtype=float)
    n = len(x)
    d = x - x.mean()
    f = np.fft.rfft(d, 2*n)
    lags = min(n - 1, 2000)
    gamma = np.fft.irfft(f*np.conj(f))[:lags]/(n - np.arange(lags))
    if gamma[0] == 0:
        return 1.0
    pairs = gamma[1:lags - 1:2] + gamma[2:lags:2]
    return n*gamma[0]/(gamma[0] + 2*(pairs*np.cumprod(pairs > 0)).sum())


def coda_ess(X):
    """effectiveSize of the R package coda, for each row of X."""
    n = X.shape[1]
    p = min(n - 1, int(np.floor(10*np.log10(n))))
    Xc = X - X.mean(axis=1, keepdims=True)
    r = np.stack([np.einsum("ij,ij->i", Xc[:, :n - k], Xc[:, k:]) for k in range(p + 1)], axis=1)/n
    a = np.zeros((len(X), 0))
    v, s = [r[:, 0]], [np.zeros(len(X))]
    for m in range(1, p + 1):
        phi = (r[:, m] - (a*r[:, m - 1:0:-1]).sum(axis=1))/v[-1]
        a = np.concatenate([a - phi[:, None]*a[:, ::-1], phi[:, None]], axis=1)
        v.append(v[-1]*(1 - phi**2))
        s.append(a.sum(axis=1))
    v, s = np.array(v).T, np.array(s).T
    order = np.argmin(n*np.log(v) + 2*np.arange(p + 1), axis=1)
    rows = np.arange(len(X))
    spec = v[rows, order]*n/(n - order - 1)/(1 - s[rows, order])**2
    return n*X.var(axis=1, ddof=1)/spec


def distance_matrices(histories, n, distances):
    """RNNI: Collienne & Gavryushkin 2021. Triplet: Brodal et al. 2013."""
    import tree_distances
    kids, words = tree_distances.decode(histories, n)

    def rf():
        _, ids = np.unique(words[:, n:2*n - 2].reshape(-1, words.shape[2]), axis=0, return_inverse=True)
        return tree_distances.rf_matrix(np.sort(ids.reshape(len(words), n - 2), axis=1).astype(np.int64))

    compute = {"rf": rf, "rnni": lambda: tree_distances.rnni_matrix(kids, words),
               "triplet": lambda: tree_distances.triplet_matrix(kids)}
    return [compute[d]() for d in distances]


def read_summary(path):
    with open(path) as f:
        positions = np.array(f.readline().split()[2:], dtype=float)
        names = f.readline().split()
        data = np.loadtxt(f, ndmin=2)
    return positions, names, data


def read_trees(path):
    with open(path, "rb") as f:
        b = f.read()
    k, n, interval, width = np.frombuffer(b, "<i8", 4).tolist()
    o = 40 + 8*k
    code = np.dtype("<u2" if width == 2 else "<u4")
    entry = np.dtype([("index", "<i4"), ("history", code, n - 1)])
    states, written = [], [[] for _ in range(k)]
    while o < len(b):
        state = int(np.frombuffer(b, "<i8", 1, o)[0])
        m = int(np.frombuffer(b, "<i4", 1, o + 8)[0])
        e = np.frombuffer(b, entry, m, o + 12)
        o += 12 + entry.itemsize*m
        for j, h in zip(e["index"].tolist(), e["history"]):
            written[j].append((state, h.tobytes()))
        states.append(state)
    return n, np.array(states), written, code


def frechet_ess(D2, rows):
    """frechetCorrelationESS, Magee et al. 2024, Eqs. 8-9."""
    n, k = len(rows), D2.shape[0]
    onehot = np.zeros((k, n))
    onehot[rows, np.arange(n)] = 1
    lagged = np.fft.irfft((np.conj(np.fft.rfft(onehot, 2*n))*np.fft.rfft(D2[:, rows], 2*n)).sum(axis=0), 2*n)[:n - 1]

    def pair_sums(r, indicator):
        before = np.cumsum(indicator, axis=1) - indicator
        return np.concatenate([[0.0], np.cumsum(2*np.einsum("lk,kl->l", D2[r], before))])

    first, last = pair_sums(rows, onehot), pair_sums(rows[::-1], onehot[:, ::-1])
    L = n - np.arange(n - 1)
    vf, vl = first[L]/(2*L*(L - 1)), last[L]/(2*L*(L - 1))
    rho = np.where(vf*vl == 0, np.nan, 0.5*(vf + vl - lagged/L)/np.sqrt(vf*vl))

    rho[0] = 1.0
    pairs = rho[0:n - 2:2] + rho[1::2]
    stop = np.flatnonzero(~(pairs > 0))
    if len(stop) and np.isnan(pairs[stop[0]]):
        return np.nan
    tau = -1 + 2*np.minimum.accumulate(pairs[:stop[0] if len(stop) else len(pairs)]).sum()
    return n/tau if tau > 0 else np.nan


def tree_ess(written, states, n, code, distances):
    """minPseudoESS, medianPseudoESS and frechetCorrelationESS, Magee et al. 2024."""
    index = {}
    tree = np.array([index.setdefault(h, len(index)) for _, h in written])
    present, rows, counts = np.unique(tree[np.searchsorted([s for s, _ in written], states, side="right") - 1],
                                      return_inverse=True, return_counts=True)
    if len(present) == 1:
        return np.full((len(distances), 3), np.nan)
    keys = list(index)
    histories = np.array([np.frombuffer(keys[t], code) for t in present], dtype=np.int64)
    out = []
    for D in distance_matrices(histories, n, distances):
        e = coda_ess(D[:, rows])
        out.append((e.min(), np.median(np.repeat(e, counts)), frechet_ess(D**2, rows)))
    return np.array(out)


def main():
    parser = argparse.ArgumentParser(description="ESS of a SINGER run.")
    parser.add_argument("-input", type=str, required=True, help="The _summary.log file of the run.")
    parser.add_argument("-output", type=str, required=True, help="Output file prefix.")
    parser.add_argument("-burnin", type=int, help="MCMC iterations to leave out. Default: -burnin of the run.")
    parser.add_argument("-tree_ess", nargs="+", default=["rf"], choices=["0", "rf", "rnni", "triplet"], help="Tree distances, or 0 for none. Default: rf.")
    args = parser.parse_args()
    np.seterr(divide="ignore", invalid="ignore")
    distances = [] if "0" in args.tree_ess else args.tree_ess
    if not args.input.endswith("_summary.log"):
        sys.exit(f"-input {args.input}: give the _summary.log file of the run.")
    prefix = args.input[:-len("_summary.log")]

    burnin = args.burnin
    if burnin is None:
        with open(f"{prefix}.ckp") as f:
            for line in f:
                w = line.split()
                if w and w[0] == "-burnin":
                    burnin = int(w[1])
    positions, names, data = read_summary(args.input)
    keep = data[:, 0] > burnin
    data = data[keep]
    column_ess = [scalar_ess(data[:, c]) for c in range(1, data.shape[1])]
    measures = ("minPseudoESS", "medianPseudoESS", "frechetCorrelationESS")
    if distances:
        import pyximport
        pyximport.install(language_level=3)
        n, states, written, code = read_trees(f"{prefix}_summary_trees.bin")
        states = states[states > burnin]
        tree = np.array([tree_ess(w, states, n, code, distances) for w in written])

    k = len(positions)
    per_position = np.array(column_ess[1:]).reshape(k, 4)
    print(f"rows after the burn-in (state > {burnin}): {len(data)}; positions: {k}")
    print(f"ESS of recombinations: {column_ess[0]:.0f}")
    print("%-34s %s" % ("quantiles over positions", "     min      25%      50%      75%      max"))
    rows = [("tbl", per_position[:, 0]), ("q25", per_position[:, 1]), ("q50", per_position[:, 2]), ("q75", per_position[:, 3])]
    rows += [(f"{m} ({d})", tree[:, c, i]) for c, d in enumerate(distances) for i, m in enumerate(measures)]
    for label, v in rows:
        print("%-34s %s" % (label, " ".join("%8.0f" % t for t in np.nanquantile(v, [0, 0.25, 0.5, 0.75, 1]))))
    unchanged = int(np.isnan(tree[:, 0, 0]).sum()) if distances else 0
    if unchanged:
        print(f"{unchanged} of the {k} positions keep one tree after the burn-in and have no tree ESS")
    with open(f"{args.output}_ess.tsv", "w") as f:
        f.write("column\tESS\n")
        for name, e in zip(names[1:], column_ess):
            f.write(f"{name}\t{e}\n")
        for c, d in enumerate(distances):
            for j in range(k):
                for i, m in enumerate(measures):
                    f.write(f"{m}_{d}_{j}\t{tree[j, c, i]}\n")


if __name__ == "__main__":
    main()
