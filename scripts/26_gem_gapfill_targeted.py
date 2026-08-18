#!/usr/bin/env python3
"""Phase 6 -- targeted cofactor gapfilling (replaces the intractable global MILP).

The first attempt handed COBRApy's gapfill() the whole ModelSEED universal pool (~40k
reactions) for each blocked target. That is a MILP whose size is set by the pool, not by the
problem, and it ran 10.5 h without converging on a single target.

This version sizes the problem to the question. For each blocked cofactor it walks BACKWARD
through the universal network from that metabolite, collecting only reactions that could
plausibly contribute within a few steps, and runs the MILP against that restricted pool. The
candidate set drops from ~40k to hundreds or a few thousand, which solves in seconds.

Reactions added carry no gene association and are tagged notes['gapfilled']='true', so
gapfilled chemistry can never be mistaken for gene-supported chemistry. A per-target time
limit prevents any repeat of the original hang.
"""
import argparse, collections, csv, pathlib, re, sys, time
import cobra
from cobra import Metabolite, Reaction
from cobra.flux_analysis import gapfill

ROOT = pathlib.Path(__file__).resolve().parents[1]
CPD_RE = re.compile(r"cpd\d+")
# NOTE: two of these ids were wrong in the first three scripts and the error propagated into
# the notes and manuscript. cpd00166 is Calomide (a cobalamin) and cpd00557 is Siroheme --
# neither is CoA or a chitin precursor. Correct ids verified against
# refs/modelseed/compounds.tsv: CoA = cpd00010, UDP-N-acetylglucosamine = cpd00037.
TARGETS = {"cpd00010": "CoA", "cpd00003": "NAD", "cpd00006": "NADP",
           "cpd00015": "FAD", "cpd00104": "biotin", "cpd00037": "UDP-GlcNAc (chitin precursor)"}

BIOMASS = {
    "cpd00035": 0.09, "cpd00051": 0.05, "cpd00132": 0.04, "cpd00041": 0.05,
    "cpd00084": 0.01, "cpd00053": 0.05, "cpd00023": 0.09, "cpd00033": 0.08,
    "cpd00119": 0.02, "cpd00322": 0.05, "cpd00107": 0.08, "cpd00039": 0.06,
    "cpd00060": 0.02, "cpd00066": 0.03, "cpd00129": 0.04, "cpd00054": 0.05,
    "cpd00161": 0.05, "cpd00065": 0.01, "cpd00069": 0.02, "cpd00156": 0.06,
    "cpd00002": 0.05, "cpd00038": 0.02, "cpd00052": 0.02, "cpd00062": 0.02,
    "cpd00115": 0.01, "cpd00241": 0.01, "cpd00356": 0.01, "cpd00357": 0.01,
    "cpd00010": 0.005, "cpd00003": 0.005, "cpd00006": 0.005, "cpd00015": 0.002,
}


def parse_stoich(s):
    out = []
    for term in (s or "").split(";"):
        p = term.split(":")
        if len(p) < 3: continue
        try: c = float(p[0])
        except ValueError: continue
        if CPD_RE.fullmatch(p[1]):
            out.append((c, p[1], int(p[2]) if p[2].isdigit() else 0))
    return out


def build_universal(rtsv, ctsv, skip_ids):
    names = {}
    with open(ctsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            names[row["id"]] = row.get("name", row["id"])
    uni = cobra.Model("universal")
    comp = {0: "c", 1: "e", 2: "p"}
    mets, rxns = {}, []
    with open(rtsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("is_obsolete") or "0") == "1" or row["id"] in skip_ids: continue
            st = parse_stoich(row.get("stoichiometry"))
            if not st: continue
            r = Reaction(row["id"])
            rev = row.get("reversibility") or "="
            r.lower_bound = -1000.0 if rev in ("=", "<") else 0.0
            r.upper_bound = 1000.0 if rev in ("=", ">") else 0.0
            if r.lower_bound == 0 and r.upper_bound == 0: r.upper_bound = 1000.0
            sm = {}
            for c, cid, k in st:
                key = f"{cid}_{comp.get(k,'c')}"
                m = mets.get(key) or Metabolite(key, name=names.get(cid, cid), compartment=comp.get(k, "c"))
                mets[key] = m
                sm[m] = sm.get(m, 0) + c
            sm = {m: v for m, v in sm.items() if v != 0}
            if sm:
                r.add_metabolites(sm)
                rxns.append(r)
    uni.add_reactions(rxns)
    return uni


def producers_index(uni):
    """metabolite id -> reactions that can yield it in some allowed direction."""
    idx = collections.defaultdict(set)
    for r in uni.reactions:
        for m, c in r.metabolites.items():
            if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
                idx[m.id].add(r.id)
    return idx


def backward_pool(target_id, uni, idx, have, depth, cap):
    """Reactions within `depth` backward steps of the target, capped in size."""
    frontier, seen_m, chosen = {target_id}, {target_id}, set()
    for _ in range(depth):
        nxt = set()
        for mid in frontier:
            for rid in idx.get(mid, ()):
                if rid in chosen: continue
                chosen.add(rid)
                if len(chosen) >= cap: return chosen
                for m, c in uni.reactions.get_by_id(rid).metabolites.items():
                    if m.id in seen_m or m.id in have: continue
                    seen_m.add(m.id); nxt.add(m.id)
        if not nxt: break
        frontier = nxt
    return chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT / "models/BOM_ss5_medium.xml"))
    ap.add_argument("--out", default=str(ROOT / "models/BOM_ss5_gapfilled.xml"))
    ap.add_argument("--reactions", default=str(ROOT / "refs/modelseed/reactions.tsv"))
    ap.add_argument("--compounds", default=str(ROOT / "refs/modelseed/compounds.tsv"))
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--cap", type=int, default=2500)
    ap.add_argument("--timeout", type=int, default=300, help="seconds per target")
    args = ap.parse_args()

    model = cobra.io.read_sbml_model(args.model)
    print(f"start: {len(model.reactions):,} reactions, {len(model.metabolites):,} metabolites", flush=True)

    uni = build_universal(args.reactions, args.compounds, {r.id for r in model.reactions})
    idx = producers_index(uni)
    have = {m.id for m in model.metabolites}
    print(f"universal pool: {len(uni.reactions):,} reactions (restricted per target below)", flush=True)

    added_all, report = set(), []
    for cid, label in TARGETS.items():
        key = f"{cid}_c"
        if key not in {m.id for m in model.metabolites}:
            model.add_metabolites([Metabolite(key, name=label, compartment="c")])
        met = model.metabolites.get_by_id(key)
        with model:
            dm = model.add_boundary(met, type="demand"); model.objective = dm
            if (model.slim_optimize() or 0) > 1e-6:
                print(f"  {label:<18} already producible", flush=True)
                report.append((label, "already producible", 0)); continue

        pool_ids = backward_pool(key, uni, idx, have, args.depth, args.cap)
        sub = cobra.Model("sub")
        sub.add_reactions([uni.reactions.get_by_id(r).copy() for r in pool_ids])
        t0 = time.time()
        try:
            with model:
                dm = model.add_boundary(met, type="demand"); model.objective = dm
                sol = gapfill(model, sub, demand_reactions=False, iterations=1)
            rxns = sol[0]
            for r in rxns:
                if r.id not in {x.id for x in model.reactions}:
                    c = r.copy(); c.gene_reaction_rule = ""
                    c.notes["gapfilled"] = "true"
                    model.add_reactions([c]); added_all.add(c.id)
            print(f"  {label:<18} pool={len(pool_ids):>5}  +{len(rxns)} rxn  "
                  f"({time.time()-t0:.1f}s): {', '.join(r.id for r in rxns[:5])}", flush=True)
            report.append((label, f"gapfilled +{len(rxns)}", len(pool_ids)))
        except Exception as e:
            print(f"  {label:<18} pool={len(pool_ids):>5}  FAILED after {time.time()-t0:.1f}s "
                  f"({type(e).__name__}: {str(e)[:70]})", flush=True)
            report.append((label, f"failed: {type(e).__name__}", len(pool_ids)))

    print(f"\ngapfilled reactions added: {len(added_all)}", flush=True)

    # re-test each target now that all additions are in
    ok = 0
    for cid, label in TARGETS.items():
        met = model.metabolites.get_by_id(f"{cid}_c")
        with model:
            dm = model.add_boundary(met, type="demand"); model.objective = dm
            v = model.slim_optimize() or 0
        state = "producible" if v > 1e-6 else "STILL BLOCKED"
        ok += v > 1e-6
        print(f"  {label:<18} {state}", flush=True)
    print(f"cofactors producible: {ok}/{len(TARGETS)}", flush=True)

    bio = Reaction("BIOMASS_fungal")
    bio.name = "Coarse fungal biomass (uncurated)"
    bio.lower_bound, bio.upper_bound = 0.0, 1000.0
    sm = {}
    for cid, coef in BIOMASS.items():
        try: sm[model.metabolites.get_by_id(f"{cid}_c")] = -coef
        except KeyError: pass
    bio.add_metabolites(sm)
    model.add_reactions([bio]); model.objective = bio
    v = model.slim_optimize() or 0
    print(f"\nbiomass precursors present: {len(sm)}/{len(BIOMASS)}", flush=True)
    print(f"biomass flux on MNM medium: {v:.4f}", flush=True)

    blocked = len(cobra.flux_analysis.find_blocked_reactions(model))
    print(f"blocked: {blocked:,}/{len(model.reactions):,} ({100*blocked/len(model.reactions):.1f}%) "
          f"-> {len(model.reactions)-blocked:,} can carry flux", flush=True)

    out = pathlib.Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(out))
    print(f"wrote {out}", flush=True)
    print("NOTE: gapfilled reactions have NO gene association; biomass objective is "
          "uncurated and coarse. Not a validated growth prediction.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
