# Python version of the SocialMapper reembedding and labelling

Place any pair of rats, real or generated, into the published SDANNCE action map (one rat's own movement) and social map (the pair). Get back per-frame map positions and fine and coarse cluster labels with the same numbering as the stored SDANNCE labels. No MATLAB needed.

These are ports of:
- `Scripts/socialFeaturesAndReembedding.m` (social features and reembedding);
- `AllCode/runReembedBigRun_HDZ.m` (action features and reembedding);
- `findWatershedRegions_v2.m` (map position to cluster).

They use the same training sets, principal components and watershed regions.

## Setup

1. Python 3.9+ and `pip install -r requirements.txt` (numpy 1.22 or newer, scipy, pandas, scikit-learn, matplotlib).
2. Download the training files from figshare (DOI 10.6084/m9.figshare.24796620) into `SocialMapperSampleData/` at the top of this repository:
   - `train1228.mat` (action map);
   - `trainingSocial_20230131.mat` (social map);
   - `sampleTracksSocial.mat` (example pair).
3. The PCA files (`vecsVals_2022_12_22.mat`, `vecsValsSOC_2022_12_28.mat`) are already at the top of this repository. The cluster tables, watershed images and speed mixtures are in `python/assets/`. To keep the data elsewhere, set `SOCIALMAPPER_ROOT` to a folder holding the PCA files and `SocialMapperSampleData/`.

## Run

```
cd python
python embed_pair.py --mat ../SocialMapperSampleData/sampleTracksSocial.mat --rat_keys p1,p2 --stop 3000 --out sample_maps.npz --csv
```

- Input: two arrays, one per rat, world coordinates in mm, 50 fps, rat23 joint order (SDANNCE m1/m2), shaped (T, 23, 3) or (T, 3, 23).
- Reading the input: `--mat` or `--npz` with `--rat_keys` names both arrays, or `--npy a.npy,b.npy`.
- `--index k` picks one trajectory from (K, T, 23, 3) arrays.
- `--start`, `--stop` and `--keep_from` select frames.
- `--no_action` or `--no_social` runs one map only.
- `--csv` also writes one row per frame.

From Python:

```
import embed_pair as ep
out = ep.embed_pair(rat1, rat2)          # dict of per-frame arrays, both maps, both rats
```

## Output (per frame, for each rat r = 1, 2)

| Name | Meaning |
| --- | --- |
| `social_fine{r}` | fine social cluster, 1-156, from rat r's point of view |
| `social_coarse{r}` | 0 no interaction; 1 rat 1 engaged; 3 rat 1 lightly engaged; 2 rat 2 engaged; 4 rat 2 lightly engaged; 7 light mutual interaction; 8 mutual interaction ("rat 1" = rat r) |
| `action_fine{r}` | fine action cluster, 1-162 |
| `action_coarse{r}` | 1 idle, 2 sniff/head, 3 groom, 4 scrunched, 5 active crouched, 6 reared, 7 explore, 8 locomotion, 9 fast/error |
| `social_z{r}`, `action_z{r}` | map positions (social map spans -75 to 75, action map -58 to 58) |
| `social_inhull{r}`, `action_inhull{r}` | False when a frame lands outside the part of the map its training neighbours cover |
| `action_nn_fine{r}` | 1-nearest-neighbour class on the same features (no t-SNE placement, no time segmentation) |

Every fine cluster's description is in `assets/social_vocab.csv` and `assets/action_vocab.csv`. The `watershed_id` column there is the region number drawn on the maps; `position` is the cluster number the code returns.

## Accuracy

Each check below runs the same frames through the original MATLAB, or compares with the stored SDANNCE positions and labels.

| Check | Social map | Action map |
| --- | --- | --- |
| Python vs MATLAB, features | max difference 6e-5 | heights and speeds 1e-14; wavelets max 0.25 (log amplitude) |
| Python vs MATLAB, map position (median) | under 0.001 units for in-hull frames | 0.009 units |
| Python vs stored SDANNCE positions (median) | 0.08-0.12 units of 150 | 0.006 units of 116, from raw s-DANNCE predictions |
| Same fine cluster as stored | 89-92% of frames | 95% (97% when the stored positions are labelled the same way) |
| Same coarse class as stored | not checked | 98% |

These checks used up to a few thousand frames from 1-4 recordings each.

**Which tracks go in matters for the action map.** The stored action positions and labels were computed from the raw s-DANNCE predictions (`save_data_AVG0.mat`). Smoothed or cleaned tracks, low-passed for example to about 7% of the raw power at 10-20 Hz, move frames on the action map: in our check by a median of about 6 map units, with about a third of frames in the same fine cluster. That is because the action features read wavelets up to 20 Hz on a log scale, where sub-millimetre jitter dominates. The social features stop at 5 Hz and are unaffected. So:
- to reproduce stored action labels, give the code the tracks those labels were made from;
- to compare two sets of tracks on the action map, process both the same way and check that their power above 5 Hz is similar.

## Notes

- **Edges.** The slowest social wavelet is 0.2 Hz, about 4 s wide. Put several seconds of recording in front of a stretch of interest (then drop it with `--keep_from`), and treat the last few seconds of any trajectory with caution.
- **Body scale (action map).** Features are divided by body length, the 95th percentile of the snout to tail-base distance over the input. For short stretches, take it from the whole recording (`--scale_from` or `--scale_from_context N`).
- **Labelling settings, matched to the stored labels.**
  - Action: the full watershed image, a minimum run of 5 frames, no median filter.
  - Social: `LL2soc`, a minimum run of 10 frames, a 2-frame median.
  - The speed mixture of the pause/movement rule was fitted once on stored map positions of real recordings and is reused for any input.
- **Cost.** About 21 ms (social) and 36 ms (action) per frame per rat on one CPU core. The code is single-threaded; run files in parallel.
