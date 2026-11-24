#!/usr/bin/env python3
"""Map coordinates -> per-frame cluster label, for the action map and the social map.

One generic port of SocialMapper/visualization/findWatershedRegions_v2.m (MotionMapper pause segmentation):
  median-filter z (n = median_n), map speed, 2-component Gaussian mixture on log10 speed, pRest = 1 - P(fast component);
  frames are assigned to the watershed region they sit in when they belong to a slow segment (pRest > 0.33, region
  unchanged, >= minrest frames, segment max pRest >= 0.67) or to a run of >= minrest frames in one region.
  Unassigned frames take the label of the nearest assigned frame in time (fillmissing 'nearest', as in sDANNCE_main.m).
The watershed id is then converted to the cluster POSITION (row of the description table), which is what the stored
SDANNCE labels hold (llac / lljc), and to the coarse class.

Settings per map:
  action  R = 58 (xx = linspace(-58, 58, 501), sDANNCE_main.m), median_n 1, minrest 5 (sDANNCE_main.m settings);
          reproduces the stored action labels from the stored coordinates on 98% of frames (coarse 99%).
  social  R = 75, median_n 2, minrest 10: reproduces the stored social labels from the stored coordinates on 93-94% of frames
          (the .m default minrest 5 gives 75-78%).
The mixture is fitted once on real stored map coordinates and reused for generated data (same speed threshold).
"""
import numpy as np
import pandas as pd
from scipy.ndimage import label as cc_label, distance_transform_edt
from sklearn.mixture import GaussianMixture

NPIX, DT = 501, 0.01


def medfilt(x, n):
    """MATLAB medfilt1(x, n) for n = 1 (identity) or 2 (mean of x(k-1), x(k), zero padded)."""
    if n <= 1:
        return np.asarray(x, float)
    assert n == 2, 'only medianLength 1 or 2 are used'
    p = np.r_[0., x]
    return 0.5 * (p[:-1] + p[1:])


def speed(z, median_n):
    s = np.stack([medfilt(z[:, 0], median_n), medfilt(z[:, 1], median_n)], 1)
    v = np.r_[0, np.linalg.norm(np.diff(s, axis=0), axis=1)] / DT
    return s, v


def fit_gmm(z_list, median_n, n=50000, seed=0):
    rng = np.random.default_rng(seed)
    lv = []
    for z in z_list:
        v = speed(z, median_n)[1]
        lv.append(np.log10(v[v > 0]))
    lv = np.concatenate(lv)
    lv = rng.choice(lv, min(n, len(lv)), replace=False)[:, None]
    return GaussianMixture(2, random_state=seed).fit(lv)


def save_gmm(g, path, **extra):
    np.savez(path, weights=g.weights_, means=g.means_, covariances=g.covariances_, **extra)


def load_gmm(path):
    P = np.load(path)
    g = GaussianMixture(2, covariance_type='full')
    g.weights_, g.means_, g.covariances_ = P['weights'], P['means'], P['covariances']
    g.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(P['covariances']))
    return g


class MapLabeler:
    """LL (501 x 501 watershed image, 0 = boundary/background), R (map half-width), wid2pos (watershed id -> position),
    pos2class (position -> coarse class, index 0 unused)."""

    def __init__(self, LL, R, wid2pos, pos2class, median_n, minrest, gmm):
        self.LL, self.R, self.median_n, self.minrest, self.gmm = np.asarray(LL, int), float(R), median_n, minrest, gmm
        self.wid2pos, self.pos2class = np.asarray(wid2pos, int), np.asarray(pos2class, int)

    def pixel(self, s):
        return np.clip(np.round((s + self.R) * NPIX / (2 * self.R)).astype(int), 1, NPIX) - 1   # row = y, col = x

    def lookup(self, z):
        """Position of the region each frame sits in, no time segmentation (0 = on a boundary)."""
        v = self.pixel(np.asarray(z, float))
        return self.wid2pos[self.LL[v[:, 1], v[:, 0]]]

    def lookup_filled(self, z):
        """As lookup, but a point on a region border takes the nearest region (used for the training frames of the 1-NN class)."""
        if not hasattr(self, '_filled'):
            _, (iy, ix) = distance_transform_edt(self.LL == 0, return_indices=True)
            self._filled = self.LL[iy, ix]
        v = self.pixel(np.asarray(z, float))
        return self.wid2pos[self._filled[v[:, 1], v[:, 0]]]

    def regions(self, z, pthr=(.67, .33)):
        """Per-frame cluster position (1-based; 0 only if nothing in the trajectory could be assigned)."""
        z = np.asarray(z, float)
        s, v = speed(z, self.median_n)
        hi = int(np.argmax(self.gmm.means_.ravel()))
        pr = np.ones(len(v))
        ok = v > 0
        pr[ok] = 1 - self.gmm.predict_proba(np.log10(v[ok])[:, None])[:, hi]
        px = self.pixel(s)
        w = self.LL[px[:, 1], px[:, 0]]
        same = np.r_[True, np.diff(w) == 0]
        out = np.zeros(len(w), int)
        lab, n = cc_label((pr > min(pthr)) & same)
        if n:
            idx = np.arange(len(w))
            for _, ii in pd.Series(idx[lab > 0]).groupby(lab[lab > 0]):
                ii = ii.values
                if len(ii) >= self.minrest and pr[ii].max() >= max(pthr):
                    vv = w[ii]; vv = vv[vv > 0]
                    if len(vv):
                        out[ii] = np.bincount(vv).argmax()
        c = np.flatnonzero(np.diff(w) != 0) + 1; st = np.r_[0, c]; ln = np.diff(np.r_[st, len(w)])
        for a, l in zip(st, ln):                                          # runs of one region >= minrest frames
            if l >= self.minrest and w[a] > 0:
                out[a:a + l] = w[a]
        if (out > 0).any():
            _, (j,) = distance_transform_edt(out == 0, return_indices=True)
            out = out[j]
        return self.wid2pos[out]

    def coarse(self, pos):
        return self.pos2class[np.asarray(pos, int)]
