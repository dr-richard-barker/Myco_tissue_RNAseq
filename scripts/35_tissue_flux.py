#!/usr/bin/env python3
"""Phase 6c -- compare the context-specific tissue models.

Three questions:
  1. which reactions are retained in which tissue (presence/absence);
  2. among reactions shared with the baseline, which carry different flux;
  3. what each tissue is predicted to exchange with its surroundings, which is the
     closest the model gets to a statement about exudate composition.

"Baseline" is the two mycelial tissues. They are the appropriate comparator biologically
(undifferentiated vegetative growth) even though their marker-level signatures were not
individually defensible: contextualisation uses network-wide expression rather than individual
marker calls. Exudophore results rest on two libraries and are reported as such.
"""
import argparse, collections, csv, pathlib, sys, warnings
warnings.filterwarnings("ignore")
import cobra
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIS = ["Exuding_mycelium", "Fuzzy_mycelium", "Exudophore", "Nodule"]
BASELINE = ["Exuding_mycelium", "Fuzzy_mycelium"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(ROOT / "models/tissue"))
    ap.add_argument("--out", default=str(ROOT / "results/tissue_metabolism"))
    args = ap.parse_args()
    d = pathlib.Path(args.dir); out = pathlib.Path(args.out); out.mkdir(parents=True, exist_ok=True)

    models, flux = {}, {}
    for t in TIS:
        f = d / f"BOM_ss5_{t}.xml"
        if not f.exists():
            print(f"  missing {f.name}"); continue
        models[t] = cobra.io.read_sbml_model(str(f))
        fp = d / f"flux_{t}.csv"
        if fp.exists():
            flux[t] = pd.read_csv(fp)
    print(f"loaded {len(models)} tissue models", flush=True)

    base = cobra.io.read_sbml_model(str(ROOT / "models/BOM_ss5_gapfilled.xml"))
    rname = {r.id: (r.name or r.id) for r in base.reactions}
    rgpr = {r.id: r.gene_reaction_rule for r in base.reactions}

    # ---------- 1. presence / absence ----------
    allr = sorted({r.id for m in models.values() for r in m.reactions})
    rows = []
    for rid in allr:
        present = {t: (rid in {x.id for x in models[t].reactions}) for t in models}
        rows.append(dict(reaction=rid, name=rname.get(rid, "")[:80],
                         **{t: int(present[t]) for t in models},
                         gpr=rgpr.get(rid, "")[:60]))
    pa = pd.DataFrame(rows)
    pa.to_csv(out / "reaction_presence.csv", index=False)
    print(f"\nunion of retained reactions: {len(allr)}")
    for t in models:
        print(f"  {t:<20} {int(pa[t].sum()):>4}")

    if "Exudophore" in models:
        base_cols = [t for t in BASELINE if t in models]
        uniq = pa[(pa["Exudophore"] == 1) & (pa[base_cols].sum(axis=1) == 0)]
        shared = pa[(pa["Exudophore"] == 1) & (pa[base_cols].sum(axis=1) > 0)]
        print(f"\nExudophore reactions absent from both baseline tissues: {len(uniq)}")
        for _, r in uniq.head(25).iterrows():
            print(f"    {r['reaction']:<12} {r['name'][:66]}")
        uniq.to_csv(out / "exudophore_unique_reactions.csv", index=False)
        print(f"shared with baseline: {len(shared)}")

    # ---------- 2. flux contrast on shared reactions ----------
    if flux:
        med = {}
        for t, df in flux.items():
            med[t] = df.median(axis=0, numeric_only=True)
        M = pd.DataFrame(med)
        M.to_csv(out / "median_flux.csv")
        if "Exudophore" in M.columns:
            bl = [c for c in BASELINE if c in M.columns]
            comp = M.dropna(subset=["Exudophore"]).copy()
            comp["baseline"] = comp[bl].mean(axis=1)
            comp = comp.dropna(subset=["baseline"])
            comp["abs_diff"] = (comp["Exudophore"].abs() - comp["baseline"].abs())
            comp["name"] = [rname.get(i, "")[:70] for i in comp.index]
            comp = comp.sort_values("abs_diff", ascending=False)
            comp.to_csv(out / "flux_contrast.csv")
            print("\nHighest median |flux| in exudophore relative to baseline:")
            for i, r in comp.head(15).iterrows():
                print(f"    {i:<12} {r['name'][:56]:<58} {r['Exudophore']:>9.2f} vs {r['baseline']:>8.2f}")

    # ---------- 3. predicted exchange (what leaves the cell) ----------
    ex_rows = []
    for t, m in models.items():
        for r in m.reactions:
            if not r.boundary:
                continue
            met = list(r.metabolites)[0]
            with m:
                m.objective = r
                try:
                    vmax = m.slim_optimize()
                except Exception:
                    vmax = float("nan")
            ex_rows.append(dict(tissue=t, reaction=r.id,
                                metabolite=met.name or met.id, max_flux=vmax))
    if ex_rows:
        ex = pd.DataFrame(ex_rows)
        ex.to_csv(out / "exchange_capacity.csv", index=False)
        piv = ex.pivot_table(index="metabolite", columns="tissue", values="max_flux")
        if "Exudophore" in piv.columns:
            bl = [c for c in BASELINE if c in piv.columns]
            piv["baseline"] = piv[bl].mean(axis=1)
            piv["exudophore_gain"] = piv["Exudophore"] - piv["baseline"]
            piv = piv.sort_values("exudophore_gain", ascending=False)
            piv.to_csv(out / "secretion_contrast.csv")
            print("\nGreatest secretion/exchange capacity in exudophore vs baseline:")
            for i, r in piv.head(12).iterrows():
                if r["exudophore_gain"] == r["exudophore_gain"] and abs(r["exudophore_gain"]) > 1e-6:
                    print(f"    {str(i)[:46]:<48} {r['Exudophore']:>9.2f} vs {r['baseline']:>8.2f}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
