# Which model file to use

Four files per reference (`BOM_ss5_*`, `PC9.15_*`), from three pipeline stages —
**not** three independent variants. `_medium.xml` is a superseded intermediate,
not an alternative to `_gapfilled.xml`.

| File | Stage | Biomass objective | Use for |
|---|---|---|---|
| `*_draft.xml` | `14_build_gem.py` — EC-mapped scaffold | none | nothing; 81.5% of reactions blocked (NOTES.md §10) |
| `*_medium.xml` | `16_gem_medium.py` — adds transport reactions, constrains exchanges to MNM v3 | **none, by design** | connectivity diagnostics only (which biomass precursors are reachable) — explicitly documented as "not a validated growth prediction" (NOTES.md §10b). This file is the *input* to the next stage, not a finished product. |
| `*_gapfilled.xml` | `26_gem_gapfill_targeted.py` — reads `*_medium.xml`, gapfills the 3 remaining blocked cofactor precursors, then adds `BIOMASS_fungal` and sets it as the objective | **yes** (`BIOMASS_fungal`, coarse/uncurated) | **the model to use for any flux prediction.** Already on the MNM v3 medium (that constraint carries over from `_medium.xml`) — do not look for a separate "final medium-constrained" file, this is it. Biomass flux = 125.0 for both references (NOTES.md §14c). |

`cobra.io.read_sbml_model()` on a `*_medium.xml` file will report "No objective
coefficients in model" and any FBA run against it is meaningless — that is
expected, not a bug. If you need the medium constraints *and* a working
objective, use `*_gapfilled.xml`.

(Added after an external consumer of this repo — the `LunarFarm-BLSS`
minimal-consortium module — briefly misread the missing objective in
`*_medium.xml` as a build defect. See that repo's `consortium/FINDINGS.md`
for the corrected story.)
