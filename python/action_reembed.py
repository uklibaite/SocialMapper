#!/usr/bin/env python3
"""Python port of the SDANNCE ACTION-map reembedding, AllCode/runReembedBigRun_HDZ.m.

That script made the stored action-map coordinates (czx, czy) and, through findWatershedRegions_v2, the stored action
labels (llac) of every SDANNCE recording (runRE_bigRunLone_20221209.sh: train1228.mat, vecsVals_2022_12_22.mat).

    p      : (T, 3, 23)  one rat, world mm, 50 fps, rat23 joint order (SDANNCE m1/m2)
    z      : (T, 2)      action-map position, map spans [-58, 58]
    inhull : (T,)        inside the convex hull of the frame's training neighbours (out-of-hull z = the guess, as in the .m)
    nn_idx : (T,)        index of the closest training frame (cityblock), as the .m's cGuesses step. Its class is the map
                         region of that training frame's map position (embed_pair.py); NOT train1228.mat's CC column,
                         which uses another cluster numbering (values up to 197; the map has 163 regions)

Features (421 = the columns of train1228.mat cD):
  joints median-filtered (3) and smoothed (3); body scale sv = 250 / 95th percentile of the snout-tail-base distance;
  506 ordered inter-joint distances x sv -> minus mus -> 15 PCs (vecs) -> Morlet wavelets 0.5-20 Hz, 25 scales,
  omega0 5 -> log, floor -5  [375];  0.1 x joint heights (re-smoothed, + floor, x sv) [23];  0.5 x joint speeds
  (medfilt 10, x sv, capped at 5) [23].
  Note: runReembedBigRun_HDZ.m ADDS the floor (p1z + floorval), individualFeaturesAndReembedding.m subtracts it; the
  big-run sign (+) made the stored coordinates and is the default (floor_sign); on real data the two signs give near-identical maps.
Reembedding: findTDistributedProjections_fmin (perplexity 32, 200 neighbours), shared with the social port
(social_reembed.py, validated against MATLAB).

sv and the floor are computed over the trajectory you pass in, as the big run did over each 30-min recording. For a
short or generated trajectory pass them explicitly (sv=, floor=) from the real recording it was seeded from.
EDGES: the slowest wavelet (0.5 Hz) spans about 2 s; the first and last ~1-2 s of any trajectory are edge-affected.
"""
import os
import numpy as np
import scipy.io as sio
from scipy.spatial.distance import cdist
import social_reembed as sr

ASSETS = os.environ.get('SOCIALMAPPER_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # repository root
TRAIN = os.path.join(ASSETS, 'SocialMapperSampleData', 'train1228.mat')
VECS = os.path.join(ASSETS, 'vecsVals_2022_12_22.mat')
MINF, MAXF, NPC, FLOOR_CLAMP = 0.5, 20.0, 15, -5.0
_A = {}


def _aux():
    if not _A:
        v = sio.loadmat(VECS)
        _A['vecs'] = v['vecs'][:, :NPC]; _A['mus'] = v['mus'].ravel()
        t = sio.loadmat(TRAIN)
        _A['cD'] = t['cD']; _A['ydata'] = t['ydata']; _A['CC'] = t['CC'].ravel().astype(int)
    return _A


def freqs(minf=MINF, maxf=MAXF, n=sr.NUM_PERIODS):
    minT, maxT = 1.0 / maxf, 1.0 / minf
    Ts = minT * 2.0 ** (np.arange(n) * np.log(maxT / minT) / (np.log(2) * (n - 1)))
    return (1.0 / Ts)[::-1]


def body_scale(p):
    """sv = 250 / p95(snout - tail base distance) and the floor (10th percentile of the two feet heights), from
    smoothed joints. MATLAB prctile = numpy 'hazen'."""
    ma = smooth_joints(p)
    sj = np.linalg.norm(ma[:, :, 0] - ma[:, :, 6], axis=1)
    return 250.0 / np.percentile(sj, 95, method='hazen'), float(np.percentile(ma[:, 2][:, [18, 22]].ravel(), 10, method='hazen'))


def smooth_joints(p):
    p = np.asarray(p, float); ma = np.empty_like(p)
    for y in range(3):
        for x in range(23):
            ma[:, y, x] = sr._ms3(p[:, y, x])
    return ma


def action_features(p, sv=None, floor=None, floor_sign=+1.0):
    """(T, 3, 23) -> (T, 421)."""
    A = _aux()
    ma = smooth_joints(p)
    T = ma.shape[0]
    if sv is None or floor is None:
        sv0, f0 = body_scale(p)
        sv = sv0 if sv is None else sv
        floor = f0 if floor is None else floor
    Xi, Yi = np.meshgrid(np.arange(23), np.arange(23))
    Xi, Yi = Xi.flatten('F'), Yi.flatten('F'); keep = Xi != Yi
    D = np.linalg.norm(ma[:, :, Xi[keep]] - ma[:, :, Yi[keep]], axis=1) * sv          # (T, 506)
    proj = (D - A['mus']) @ A['vecs']                                                  # (T, 15)
    f = freqs()
    W = np.concatenate([sr._fast_wavelet(proj[:, i], f).T for i in range(NPC)], axis=1)
    with np.errstate(divide='ignore'):
        W = np.log(W)
    W[W < FLOOR_CLAMP] = FLOOR_CLAMP
    jv = np.zeros((T, 23))
    for i in range(23):
        step = np.sqrt(np.sum(np.diff(ma[:, :, i], axis=0) ** 2, axis=1))
        jv[:, i] = np.r_[0.0, sr._medfilt1(step, 10)]
    jv = jv * sv; jv[jv >= 5] = 5.0
    pz = np.stack([sr._ms3(ma[:, 2, i]) for i in range(23)], 1)
    pz = (pz + floor_sign * floor) * sv
    return np.concatenate([W, 0.1 * pz, 0.5 * jv], axis=1)


def reembed(X, batch=500, verbose=False):
    """findTDistributedProjections_fmin against train1228 -> z (T, 2), inhull (T,); plus nn_idx (T,), the closest training frame."""
    A = _aux()
    sd, yd = A['cD'], A['ydata']
    N = len(X)
    z = np.zeros((N, 2)); guess = np.zeros((N, 2)); inch = np.zeros(N, bool); nn = np.zeros(N, int)
    opt = dict(method='Nelder-Mead', options=dict(maxiter=100, xatol=1e-4, fatol=1e-4, disp=False))
    from scipy.optimize import minimize
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path as MplPath
    for b0 in range(0, N, batch):
        b1 = min(b0 + batch, N)
        D = cdist(X[b0:b1], sd)                                                     # euclidean, as pdist2
        nn[b0:b1] = np.argmin(cdist(X[b0:b1], sd, 'cityblock'), axis=1)
        for k in range(b1 - b0):
            i = b0 + k
            order, p = sr._sigma_p(D[k])
            zt = yd[order]
            g = np.array([(p[:, None] * zt).sum(0), zt[int(np.argmax(p))]])
            sols = np.zeros((2, 2)); cost = np.zeros(2)
            for j in range(2):
                r = minimize(sr._kl_cost, g[j], args=(zt, p), **opt); sols[j] = r.x; cost[j] = r.fun
            try:
                inside = MplPath(zt[ConvexHull(zt).vertices]).contains_points(sols)
            except Exception:
                inside = np.array([False, False])
            if inside.any():
                cand = np.where(inside)[0]; mI = cand[int(np.argmin(cost[cand]))]; inch[i] = True
            else:
                mI = int(np.argmin(cost))
            z[i] = sols[mI]; guess[i] = g[mI]
        if verbose:
            print(f'    action reembed {b1}/{N}', flush=True)
    z[~inch] = guess[~inch]
    return z, inch, nn
