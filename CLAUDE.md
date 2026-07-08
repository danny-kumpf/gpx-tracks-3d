# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

Takes a GPX track (a recorded hike/ride) and an STL terrain model of the same area, and carves the track into the surface of the model as a groove, producing a new STL ready for 3D printing. This replaced a manual Fusion 360 workflow, which is still documented in the comment block at the top of `main.py`.

## Running it

Dependencies are managed by `uv` (`pyproject.toml` + `uv.lock`). There is no test or lint setup and no CI — this is a set of scripts, not a package.

```
uv sync              # install
uv run main.py       # the whole pipeline
uv run fix_mesh.py   # standalone pymeshfix repair, hardcoded paths
```

`requires-python = ">=3.10"`, and the code runs fine on 3.14, though `.idea/misc.xml` still pins the IDE to 3.10. Note that the `numpy-stl` dependency imports as `stl`.

`main.py` has **no CLI arguments**. To change inputs, edit the `json_inputs = load_from_json(...)` path at `main.py:77` to point at a different `data/<name>/config.json`. Cut geometry (`cut_radius_mm`, `dist_to_refine`, `n_sub`) is likewise edited in place just below it.

**A run requires network access and blocks on GUI windows.** `get_model_origin_lla` queries the OpenTopoData API for the ground elevation at the model's lower-left corner, and raises if it fails. Two PyVista windows open per GPX track and the pipeline halts until each is closed — this is by design, not a hang (see "Interactive alignment" below).

## config.json

Each `data/<name>/` holds the GPX, the source STL, and a `config.json`:

```json
{
  "in_stl_path": "...", "out_stl_path": "...",
  "gpx_path": "..."  or  ["...", "..."],
  "box_upper_right": [lon, lat, alt],
  "box_lower_left":  [lon, lat, alt],
  "hike_type": "OUT_AND_BACK" | "LOOP" | "NONE"
}
```

The two box corners are the real-world lat/lon bounds that the STL tile covers — they are how the code recovers the model's scale and georeference, so a wrong box shifts or stretches the track. Paths may be absolute or relative to the repo root; existing configs mix both. `gpx_path` accepts a single string or a list; `main.py` normalizes to a list, and multiple tracks are cut into the same output mesh.

`hike_type` drives `clean_up_track` in `main.py`: `OUT_AND_BACK` truncates the track after its highest-elevation point (so the return leg doesn't double-cut the outbound one), `LOOP` appends a copy of the first point to close the loop, `NONE` does neither.

## Coordinate pipeline

The core of the project, spread across `gpx_coordinate_transforms.py`. Four frames, in order:

1. **LLA** — straight from the GPX. Ordered **(lon, lat, alt)**, not (lat, lon). This ordering is used throughout, including the config's box corners; `highest_point_lat_lon` returns them swapped back to (lat, lon) for printing only.
2. **ECEF** — via `pyproj` (`epsg:4979` → `epsg:4978`).
3. **ENU** — meters, east/north/up, relative to a computed origin.
4. **Model mm** — called "mmm" in the code. ENU scaled by `mmm_per_m`.

Two quantities tie the frames together, both derived rather than configured:

- `get_model_scale` compares the diagonal of the real-world box (meters) against the diagonal of the STL's XY bounding box (mm) to get `mmm_per_m`.
- `get_model_origin_lla` places model `[0,0,0]` at the STL's lower-left corner. It fetches that corner's real ground elevation from OpenTopoData, then subtracts the STL's own local height there (`highest_z_near_origin`) so the origin sits at the model's base plane rather than its surface.

Elevation from these two sources is only approximate, so Z is not trusted: `move_closer_to_surface` discards each track point's computed Z entirely and snaps it to the Z of the mesh vertex nearest in **XY only**. Errors in the origin elevation therefore wash out, but errors in the box corners (which affect X/Y) do not.

## Interactive alignment

`plot_mesh_with_track_w_sliders` in `plotting.py` opens a PyVista window mid-transform with Shift X/Y/Z and uniform Scale sliders, letting you nudge the track onto the terrain by eye when the derived transform is off. On window close it applies the final shift and scale by **mutating the caller's `track_points` array in place** (`track_points[:] = ...`). There is no persistence: adjustments are lost when the process exits.

The Scale slider exists because `get_model_scale` is systematically wrong, which is the subject of its `TODO: suspicious that this doesn't agree with the logfile`. The cause: TouchTerrain rounds the requested lat/lon region outward to whole DEM cells before projecting to UTM, so the STL spans slightly *more* ground than the config's box corners describe. `get_model_scale` measures the box diagonal, underestimates the true ground extent, and returns an `mmm_per_m` that is correspondingly too large — so the track is drawn a few percent too big. Verified against the TouchTerrain logfiles: their `map scale` line equals `(cells × cell_size) diagonal / model diagonal` exactly, for both `south_arapaho` (140.688 m/mm, code says 137.399, +2.39%) and `mr_dlp_ironman` (207.976 m/mm, code says 206.991, +0.48%). The error scales with how much rounding the cell size forces, so it is worst on small models.

A real fix would derive `mmm_per_m` from the logfile's `map scale`, or reproject the corners to the tile's UTM zone rather than using an ENU chord. Until then the slider is the workaround, and the correct setting is `script_scale / logfile_scale`.

`plot_mesh_with_track` is the plain non-interactive viewer, called from `main.py` after cleanup.

## Cutting

`auto_cut.py` does the geometry work in two passes over the mesh, both driven by `cKDTree`s built from the densified track:

1. `upsample_along_track` recursively subdivides (1 triangle → 4, `n_sub` times) only the triangles with a vertex within `dist_to_refine` of a track, returning the refined mesh plus a boolean mask of which triangles were subdivided.
2. `cut_along_track` walks that mask and lowers each refined vertex's Z by `polynomial_drop(distance_to_track, cut_radius_mm)` — a smoothed half-circle profile that flattens to zero slope at the groove edge, which prints better than the sharp `half_circle_drop`. Untouched triangles are skipped because the mask already proved them too far away.

Both take a **list** of KD-trees, one per GPX track, and take the minimum distance across all of them.

Note that `auto_cut.py`'s own `__main__` block is stale: it passes a single `cKDTree` where the functions now expect a list, so running `python auto_cut.py` directly raises. Use it as a reference for the standalone CSV-driven flow, not as a working entry point.

## Repo conventions

`data/` and `old_data/` hold large binary STL/GPX inputs and are mostly untracked — only the five `.py` files are in git. Generated `*_autocut.stl` outputs live alongside their inputs. When adding a new model, create `data/<name>/` with the STL, GPX, and a `config.json`, then point `main.py` at it.
