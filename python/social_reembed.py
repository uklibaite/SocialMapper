#!/usr/bin/env python3
"""Faithful Python port of SocialMapper/socialFeaturesAndReembedding.m.

Turns a pair of rat tracks (p1, p2) into the 2D social-map trajectory z1/z2 by
(1) building the 194-dim social feature matrix and (2) reembedding it into the
pretrained joint social t-SNE (findTDistributedProjections_fmin). No MATLAB.

    p1, p2 : (T, 3, 23)   world-frame mm, joint order = SDANNCE m1/m2 (rat23)
    z1     : (T, 2)       animal-1 perspective  (ma1=p1, ma2=p2)  [-58, 58] space
    z2     : (T, 2)       animal-2 perspective  (ma1=p2, ma2=p1)

Feature layout (194) matches trainingSocial_20230131.mat exactly:
    [150 wavelet | 7 rat1z*.1 | 7 rat2z*.1 | 6 p1v | 6 p2v | 6 ang*10 | 6 PJ*.05 | 6 PJN*2]

EDGE BUFFER: the wavelets run down to minF=0.2 Hz (Morlet support ~4 s each side)
and velocity is medfilt over 10 frames, so the first/last ~1-2 s of any trajectory
are edge-corrupted. Reembed whole, then trim/fade the tail for display.

Aux (SocialMapper): vecsValsSOC_2022_12_28.mat (vecs, mus),
                    SocialMapperSampleData/trainingSocial_20230131.mat (trainingData, trainingEmbedding)
"""
from __future__ import annotations

import os
import numpy as np
import scipy.io as sio
from scipy.optimize import minimize
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from matplotlib.path import Path as MplPath

SM = os.environ.get('SOCIALMAPPER_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # repository root: vecsValsSOC_2022_12_28.mat here, training data in SocialMapperSampleData/

# ---- run parameters (socialFeaturesAndReembedding.m + setRunParameters defaults) ----
FS = 50.0
MINF, MAXF = 0.2, 5.0
NUM_PERIODS = 25
OMEGA0 = 5.0
PERPLEXITY = 32.0
MAX_NEIGHBORS = 200
SIGMA_TOL = 1e-5

# joint indices (0-based) — useID = [1 4 5 6 7 8 11 12 15 16 19 20 23] (MATLAB 1-based)
USEID = np.array([1, 4, 5, 6, 7, 8, 11, 12, 15, 16, 19, 20, 23]) - 1


# ------------------------------------------------------------------ small filters
def _medfilt1(x, n):
    """MATLAB medfilt1(x,n): zero-padded, even n averages the two middle values."""
    x = np.asarray(x, float)
    N = len(x)
    hr = (n - 1) // 2                       # right pad; MATLAB even-n window is left-heavy
    hl = n - 1 - hr                         # left pad = n/2 (even), (n-1)/2 (odd)
    xp = np.concatenate([np.zeros(hl), x, np.zeros(hr)])
    return np.array([np.median(xp[i:i + n]) for i in range(N)])


def _smooth3(x):
    """MATLAB smooth(x,3): 3-pt moving average, endpoints unchanged."""
    x = np.asarray(x, float)
    y = x.copy()
    if len(x) >= 3:
        y[1:-1] = (x[:-2] + x[1:-1] + x[2:]) / 3.0
    return y


def _ms3(x):
    return _smooth3(_medfilt1(x, 3))


# ------------------------------------------------------------------ Morlet wavelet
def _morlet_freqs():
    minT, maxT = 1.0 / MAXF, 1.0 / MINF
    Ts = minT * 2.0 ** (np.arange(NUM_PERIODS) * np.log(maxT / minT) / (np.log(2) * (NUM_PERIODS - 1)))
    return (1.0 / Ts)[::-1]          # fliplr(1./Ts)


def _fast_wavelet(x, f, omega0=OMEGA0, dt=1.0 / FS):
    """Port of fastWavelet_morlet_convolution_parallel.m -> amp (L, N0)."""
    x = np.asarray(x, float).ravel()
    N0 = len(x)
    odd = (N0 % 2 == 1)
    if odd:
        x = np.concatenate([x, [0.0]])
    N = len(x)                        # even
    M = N
    x = np.concatenate([np.zeros(N // 2), x, np.zeros(N // 2)])
    N = len(x)
    L = len(f)
    scales = (omega0 + np.sqrt(2 + omega0 ** 2)) / (4 * np.pi * f)
    Omegavals = 2 * np.pi * np.arange(-N // 2, N // 2) / (N * dt)
    xHat = np.fft.fftshift(np.fft.fft(x))
    if odd:                           # MATLAB idx (1-based) -> 0-based slice
        i0, i1 = M // 2, M // 2 + M - 1          # (M/2+1):(M/2+M-1)
    else:
        i0, i1 = M // 2, M // 2 + M              # (M/2+1):(M/2+M)
    amp = np.zeros((L, i1 - i0))
    corr = np.pi ** -0.25 * np.exp(0.25 * (omega0 - np.sqrt(omega0 ** 2 + 2)) ** 2)
    for i in range(L):
        m = np.pi ** -0.25 * np.exp(-0.5 * (-Omegavals * scales[i] - omega0) ** 2)
        q = np.fft.ifft(m * xHat) * np.sqrt(scales[i])
        q = q[i0:i1]
        amp[i, :] = np.abs(q) * corr / np.sqrt(2 * scales[i])
    return amp


def _find_wavelets(proj):
    """findWavelets(proj, 9, params): mode-major (mode0's 25 periods, mode1's, ...)."""
    f = _morlet_freqs()
    T, nmodes = proj.shape
    out = np.zeros((T, nmodes * NUM_PERIODS))
    for i in range(nmodes):
        out[:, i * NUM_PERIODS:(i + 1) * NUM_PERIODS] = _fast_wavelet(proj[:, i], f).T
    return out


def _ang_diff(vec1, vec2):
    """findAngleDiff: acos of the normalized dot, row-wise."""
    num = np.sum(vec1 * vec2, axis=1)
    den = np.linalg.norm(vec1, axis=1) * np.linalg.norm(vec2, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = num / den
    return np.real(np.arccos(np.clip(c, -1.0, 1.0)))


# ------------------------------------------------------------------ features
_AUX = {}


def _load_aux():
    if _AUX:
        return _AUX
    vv = sio.loadmat(os.path.join(SM, "vecsValsSOC_2022_12_28.mat"))
    _AUX["vecs6"] = vv["vecs"][:, :6]
    _AUX["mus"] = vv["mus"].ravel()
    tr = sio.loadmat(os.path.join(SM, "SocialMapperSampleData", "trainingSocial_20230131.mat"))
    _AUX["trainingData"] = tr["trainingData"]
    _AUX["trainingEmbedding"] = tr["trainingEmbedding"]
    return _AUX


# pairwise-distance index set (MATLAB meshgrid + colon flatten, exclude diagonal)
def _pair_idx():
    xIdx = np.concatenate([USEID, USEID + 23])          # 26 indices into the 46-joint stack
    Xi, Yi = np.meshgrid(xIdx, xIdx)                    # indexing='xy' == MATLAB meshgrid
    Xi, Yi = Xi.flatten("F"), Yi.flatten("F")          # (:) column-major
    keep = Xi != Yi
    return Xi, Yi, keep


def social_features(p1, p2):
    """Animal-1-perspective 194-dim feature matrix (ma1=p1, ma2=p2). p1,p2:(T,3,23)."""
    A = _load_aux()
    p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float)
    T = p1.shape[0]
    # minimal median filter + smooth per (coord, joint)
    ma1 = np.empty_like(p1)
    ma2 = np.empty_like(p2)
    for y in range(3):
        for x in range(23):
            ma1[:, y, x] = _ms3(p1[:, y, x])
            ma2[:, y, x] = _ms3(p2[:, y, x])
    maS = np.concatenate([ma1, ma2], axis=2)           # (T,3,46)

    allz = maS[:, 2][:, [18, 22, 41, 45]]              # z of joints [19 23 42 46] (1-based)
    fz = np.percentile(allz.ravel(), 10)

    # inter-keypoint distances -> PCA(6) -> wavelets
    Xi, Yi, keep = _pair_idx()
    P = np.empty((len(Xi), T))
    for i in range(len(Xi)):
        d = np.linalg.norm(maS[:, :, Xi[i]] - maS[:, :, Yi[i]], axis=1)
        P[i] = _ms3(d)
    p1Dist = P[keep].T                                 # (T, 650)
    proj = (p1Dist - A["mus"]) @ A["vecs6"]            # (T, 6)
    data2 = np.log(_find_wavelets(proj))
    data2[data2 < -3] = -3.0

    # velocities (per joint, medfilt1(.,10)), clipped at 5
    def vel(m):
        v = np.zeros((T, 23))
        for i in range(23):
            step = np.sqrt(np.sum(np.diff(m[:, :, i], axis=0) ** 2, axis=1))
            v[:, i] = np.concatenate([[0.0], _medfilt1(step, 10)])
        v[v > 5] = 5.0
        return v
    p1v, p2v = vel(ma1), vel(ma2)

    n1, s1, c1, t1 = ma1[:, :, 0], ma1[:, :, 3], ma1[:, :, 4], ma1[:, :, 6]
    n2, s2, c2, t2 = ma2[:, :, 0], ma2[:, :, 3], ma2[:, :, 4], ma2[:, :, 6]
    distc = np.linalg.norm(c1 - c2, axis=1)

    ang = np.zeros((T, 6))
    ang[:, 0] = _ang_diff(n2 - s1, n1 - s1); ang[:, 1] = _ang_diff(c2 - s1, n1 - s1)
    ang[:, 2] = _ang_diff(t2 - s1, n1 - s1); ang[:, 3] = _ang_diff(n1 - s2, n2 - s2)
    ang[:, 4] = _ang_diff(c1 - s2, n2 - s2); ang[:, 5] = _ang_diff(t1 - s2, n2 - s2)
    ang = np.pi - ang

    d2 = np.stack([np.linalg.norm(n1 - n2, axis=1), np.linalg.norm(n1 - c2, axis=1),
                   np.linalg.norm(n1 - t2, axis=1)], axis=1)          # nose1 -> rat2
    d1 = np.stack([np.linalg.norm(n2 - n1, axis=1), np.linalg.norm(n2 - c1, axis=1),
                   np.linalg.norm(n2 - t1, axis=1)], axis=1)          # nose2 -> rat1
    PJ = np.concatenate([d1, d2], axis=1)
    PJN = np.concatenate([d1 / distc[:, None], d2 / distc[:, None]], axis=1)

    rat1z = ma1[:, 2][:, [0, 3, 4, 7, 10, 11, 14]] - fz               # joints [1 4 5 8 11 12 15]
    rat2z = ma2[:, 2][:, [0, 3, 4, 7, 10, 11, 14]] - fz
    vsel = [0, 10, 11, 14, 18, 22]                                    # joints [1 11 12 15 19 23]

    nn = np.concatenate([data2, 0.1 * rat1z, 0.1 * rat2z,
                         p1v[:, vsel], p2v[:, vsel], 10 * ang,
                         0.05 * PJ, 2 * PJN], axis=1)
    return nn


# ------------------------------------------------------------------ reembedding
def _sigma_p(ds, perplexity=PERPLEXITY, tol=SIGMA_TOL, max_neighbors=MAX_NEIGHBORS):
    """returnCorrectSigma_sparse: order of the k nearest + their perplexity-matched p."""
    high = ds.max(); low = 1e-10
    order = np.argsort(ds)[:max_neighbors]
    dsn = ds[order]
    d2 = dsn ** 2

    def perp(sig):
        p = np.exp(-0.5 * d2 / sig ** 2)
        s = p.sum()
        if s <= 0:
            return 0.0, p
        p = p / s
        nz = p > 0
        H = -np.sum(p[nz] * np.log2(p[nz]))
        return 2 ** H, p

    sig = 0.5 * (high + low)
    P, p = perp(sig)
    it = 0
    while abs(P - perplexity) >= tol and it < 10000:
        if P > perplexity:
            high = sig
        else:
            low = sig
        sig = 0.5 * (high + low)
        P, p = perp(sig)
        it += 1
    return order, p


def _kl_cost(x, z, p):
    d = ((z - x) ** 2).sum(1)
    return np.log(np.sum(1.0 / (1.0 + d))) + np.sum(p * np.log(1.0 + d))


def reembed(dataMat, verbose=False, batch=1000):
    """findTDistributedProjections_fmin -> z (N,2), inConvHull (N,)."""
    A = _load_aux()
    sig_data = A["trainingData"]
    yData = A["trainingEmbedding"]
    N = dataMat.shape[0]
    z = np.zeros((N, 2))
    guesses = np.zeros((N, 2))
    inCH = np.zeros(N, bool)
    opt = dict(method="Nelder-Mead",
               options=dict(maxiter=100, xatol=1e-4, fatol=1e-4, disp=False))
    for b0 in range(0, N, batch):
        b1 = min(b0 + batch, N)
        D = cdist(dataMat[b0:b1], sig_data)              # (nb, Ntrain) euclidean == pdist2
        for k in range(b1 - b0):
            i = b0 + k
            order, p = _sigma_p(D[k])
            zt = yData[order]
            a = (p[:, None] * zt).sum(0)                    # weighted mean
            gmax = zt[int(np.argmax(p))]                    # peak neighbor
            g = np.array([a, gmax])
            sols = np.zeros((2, 2)); cost = np.zeros(2)
            for j in range(2):
                r = minimize(_kl_cost, g[j], args=(zt, p), **opt)
                sols[j] = r.x; cost[j] = r.fun
            try:
                hull = ConvexHull(zt)
                poly = MplPath(zt[hull.vertices])
                inside = poly.contains_points(sols)
            except Exception:
                inside = np.array([False, False])
            if inside.any():
                cand = np.where(inside)[0]
                mI = cand[int(np.argmin(cost[cand]))]
                inCH[i] = True
            else:
                mI = int(np.argmin(cost))
            z[i] = sols[mI]; guesses[i] = g[mI]
        if verbose:
            print(f"    reembed {b1}/{N}")
    z_out = z.copy()
    z_out[~inCH] = guesses[~inCH]                           # out-of-hull -> guess (matches .m)
    return z_out, inCH


def social_z1(p1, p2, verbose=False):
    """Animal-1 perspective social-map trajectory for a pair. p1,p2:(T,3,23) -> (T,2)."""
    return reembed(social_features(p1, p2), verbose=verbose)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self_test", type=int, default=0, help="reembed first N frames of the sample recording")
    ap.add_argument("--ref", default=None, help="MATLAB ref .mat (dataMat1, z1) to compare against")
    args = ap.parse_args()
    if args.self_test:
        s = sio.loadmat(os.path.join(SM, "SocialMapperSampleData", "sampleTracksSocial.mat"))
        n = args.self_test
        p1, p2 = s["p1"][:n], s["p2"][:n]
        F = social_features(p1, p2)
        print(f"features {F.shape}  finite={np.isfinite(F).all()}")
        z, inCH = reembed(F, verbose=True)
        print(f"z1 {z.shape}  inConvHull {inCH.mean()*100:.0f}%  "
              f"x[{z[:,0].min():.1f},{z[:,0].max():.1f}] y[{z[:,1].min():.1f},{z[:,1].max():.1f}]")
        if args.ref and os.path.exists(args.ref):
            r = sio.loadmat(args.ref)
            dm = r["dataMat1"][:n]; zr = r["z1"][:n]
            print(f"[feat] max|py-mat| = {np.abs(F-dm).max():.3e}   mean = {np.abs(F-dm).mean():.3e}")
            print(f"[z1]   median|py-mat| = {np.median(np.linalg.norm(z-zr,axis=1)):.3f}  "
                  f"max = {np.linalg.norm(z-zr,axis=1).max():.3f}")
