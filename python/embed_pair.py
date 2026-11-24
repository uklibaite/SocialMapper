#!/usr/bin/env python3
"""Put a pair of rats (real or generated, e.g. model rollouts) on the SDANNCE action map and social map, with fine and coarse labels.
Python port of the SocialMapper reembedding and labelling (no MATLAB); same features, training sets and watershed regions.

Input: two arrays, one per rat, world coordinates in mm, 50 fps, rat23 joint order (SDANNCE m1/m2), shaped (T, 23, 3) or (T, 3, 23).

Per rat r = 1, 2 (rat r's own view):
  action_z{r} (T, 2)       action-map position (map spans +-58); action_inhull{r}
  action_fine{r} (T,)      fine action cluster = position 1..162 (the numbering of the stored SDANNCE action labels)
  action_coarse{r} (T,)    1 idle, 2 sniff/head, 3 groom, 4 scrunched, 5 active crouched, 6 reared, 7 explore, 8 locomotion, 9 fast/error
  action_nn_fine{r} (T,)   1-nearest-neighbour class on the same features (no t-SNE placement, no time segmentation)
  social_z{r} (T, 2)       social-map position, rat r as "rat 1" (map spans +-75); social_inhull{r}
  social_fine{r} (T,)      fine social cluster = position 1..156 (the numbering of the stored SDANNCE social labels)
  social_coarse{r} (T,)    0 no interaction, 1 rat 1 engaged, 3 rat 1 lightly engaged, 2 rat 2 engaged, 4 rat 2 lightly engaged,
                           7 light mutual interaction, 8 mutual interaction
Descriptions of every fine cluster: assets/action_vocab.csv, assets/social_vocab.csv. See README.md (accuracy, tracks, edges).

RUN:  python embed_pair.py --mat ../SocialMapperSampleData/sampleTracksSocial.mat --rat_keys p1,p2 --stop 3000 --out sample_maps.npz --csv
"""
import os, sys, argparse, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import social_reembed as sr
import action_reembed as ar
import map_regions as mr

AS = HERE / 'assets'
_LAB = {}


def labelers():
    if not _LAB:
        Va, Vs = pd.read_csv(AS / 'action_vocab.csv'), pd.read_csv(AS / 'social_vocab.csv')
        A, S = np.load(AS / 'action_map.npz'), np.load(AS / 'social_map.npz')
        for name, V, LL, R, med, minrest, gmm in (('action', Va, A['LL'], 58.0, 1, 5, 'gmm_action_real.npz'),
                                                  ('social', Vs, S['LL2soc'], 75.0, 2, 10, 'gmm_social_real.npz')):
            w2p = np.zeros(int(LL.max()) + 1, int); w2p[V.watershed_id.values] = V.position.values
            p2c = np.zeros(len(V) + 1, int); p2c[V.position.values] = V.coarse_class.values
            _LAB[name] = mr.MapLabeler(LL, R, w2p, p2c, med, minrest, mr.load_gmm(AS / gmm))
    return _LAB


def as_tc23(x):
    x = np.asarray(x, float)
    return x if x.shape[1] == 3 else x.transpose(0, 2, 1)


def embed_pair(p1, p2, action=True, social=True, scale=None, verbose=False):
    """p1, p2: (T, 3, 23) or (T, 23, 3) mm. scale: None or ((sv1, floor1), (sv2, floor2)). Returns dict (see module doc)."""
    p = [as_tc23(p1), as_tc23(p2)]
    lab = labelers(); out = {}
    for r in (0, 1):
        k = str(r + 1)
        if action:
            sv, fl = scale[r] if scale is not None else ar.body_scale(p[r])
            z, inch, nn = ar.reembed(ar.action_features(p[r], sv=sv, floor=fl), verbose=verbose)
            L = lab['action']; fine = L.regions(z)
            nnpos = L.lookup_filled(ar._aux()['ydata'][nn])                   # map region of the closest training frame
            out.update({f'action_z{k}': z.astype(np.float32), f'action_inhull{k}': inch, f'action_fine{k}': fine.astype(np.int16),
                        f'action_coarse{k}': L.coarse(fine).astype(np.int8), f'action_nn_fine{k}': nnpos.astype(np.int16),
                        f'action_scale{k}': np.array([sv, fl])})
        if social:
            z, inch = sr.reembed(sr.social_features(p[r], p[1 - r]))
            L = lab['social']; fine = L.regions(z)
            out.update({f'social_z{k}': z.astype(np.float32), f'social_inhull{k}': inch, f'social_fine{k}': fine.astype(np.int16),
                        f'social_coarse{k}': L.coarse(fine).astype(np.int8)})
    return out


def to_frame_table(out, frames):
    cols = {'frame': frames}
    for k, v in out.items():
        if v.ndim == 1 and len(v) == len(frames):
            cols[k] = v
        elif v.ndim == 2 and len(v) == len(frames):
            cols[k + 'x'], cols[k + 'y'] = v[:, 0], v[:, 1]
    return pd.DataFrame(cols)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--npz', help='.npz file holding both rats')
    ap.add_argument('--mat', help='.mat file holding both rats (e.g. SocialMapperSampleData/sampleTracksSocial.mat, p1 and p2)')
    ap.add_argument('--rat_keys', default='gen_a1,gen_a2', help='array names of rat 1 and rat 2 in --npz')
    ap.add_argument('--index', type=int, default=None, help='pick one rollout from (K, T, 23, 3) arrays, e.g. roll_a1/roll_a2')
    ap.add_argument('--npy', default='', help='alternative to --npz: two .npy files, comma separated (rat 1, rat 2)')
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--stop', type=int, default=None)
    ap.add_argument('--keep_from', type=int, default=0, help='drop output frames before this (relative to --start), e.g. real context')
    ap.add_argument('--scale_from', default='', help='real skeleton .npy per rat, comma separated, for the action body scale')
    ap.add_argument('--scale_from_context', type=int, default=0, help='body scale from the first N frames of the input')
    ap.add_argument('--no_action', action='store_true'); ap.add_argument('--no_social', action='store_true')
    ap.add_argument('--out', required=True); ap.add_argument('--csv', action='store_true', help='also write a per-frame CSV')
    a = ap.parse_args()
    if a.mat:
        import scipy.io as sio
        k1, k2 = a.rat_keys.split(','); Z = sio.loadmat(a.mat, variable_names=[k1, k2])
        x1, x2 = Z[k1], Z[k2]
    elif a.npz:
        Z = np.load(a.npz, allow_pickle=True); k1, k2 = a.rat_keys.split(',')
        x1, x2 = Z[k1], Z[k2]
        if a.index is not None:
            x1, x2 = x1[a.index], x2[a.index]
    else:
        f1, f2 = a.npy.split(','); x1, x2 = np.load(f1), np.load(f2)
    full = [as_tc23(x1), as_tc23(x2)]
    p1, p2 = full[0][a.start:a.stop], full[1][a.start:a.stop]
    scale = None
    if a.scale_from:
        scale = [ar.body_scale(as_tc23(np.load(f))) for f in a.scale_from.split(',')]
    elif a.scale_from_context:
        scale = [ar.body_scale(f[:a.scale_from_context]) for f in full]
    t0 = time.time()
    out = embed_pair(p1, p2, action=not a.no_action, social=not a.no_social, scale=scale)
    frames = np.arange(a.start, a.start + len(p1))
    out = {k: (v[a.keep_from:] if v.ndim >= 1 and len(v) == len(frames) else v) for k, v in out.items()}
    frames = frames[a.keep_from:]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, frames=frames, source=str(a.mat or a.npz or a.npy), **out)
    if a.csv:
        to_frame_table(out, frames).to_csv(Path(a.out).with_suffix('.csv'), index=False)
    msg = [f'{len(frames)} frames, {(time.time() - t0) / 60:.1f} min']
    for k in ('action_inhull1', 'social_inhull1'):
        if k in out:
            msg.append(f'{k} {out[k].mean():.2f}')
    print('[saved]', a.out, ' '.join(msg), flush=True)


if __name__ == '__main__':
    main()
