#!/usr/bin/env python3
"""Phase 6 -- gapfill the blocked cofactor pathways and add a biomass objective.

Six biomass precursors were unreachable in both the PC9.15 and BOM_ss5 medium-constrained
models: CoA, NAD, NADP, FAD, biotin and the chitin precursor. That is a genuine limitation of
EC-driven reconstruction -- cofactor biosynthesis steps are often annotated without a
complete EC, so the EC->reaction mapping misses them -- not a reference artefact.

Approach: build a universal reaction pool from the full ModelSEED biochemistry (not just the
EC-matched subset), then use COBRApy's gapfill to find the minimal set of reactions that
makes each blocked precursor producible on the MNM medium. Gapfilled reactions are added
WITHOUT gene associations and tagged in their notes, so they are never mistaken for
gene-supported chemistry.

A fungal biomass objective is then assembled from the precursors that are reachable.
"""
import argparse, csv, pathlib, re, sys
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

# Coarse fungal biomass: 20 aa + 4 NTP + 4 dNTP + cofactors, unit-ish coefficients.
# Deliberately not a fitted composition -- it exists so the model has an objective and can be
# tested for growth, not to make quantitative yield predictions.
BIOMASS = {
    "cpd00035": 0.09, "cpd00051": 0.05, "cpd00132": 0.04, "cpd00041": 0.05,
    "cpd00084": 0.01, "cpd00053": 0.05, "cpd00023": 0.09, "cpd00033": 0.08,
    "cpd00119": 0.02, "cpd00322": 0.05, "cpd00107": 0.08, "cpd00039": 0.06,
    "cpd00060": 0.02, "cpd00066": 0.03, "cpd00129": 0.04, "cpd00054": 0.05,
    "cpd00161": 0.05, "cpd00065": 0.01, "cpd00069": 0.02, "cpd00156": 0.06,
    "cpd00002": 0.05, "cpd00038": 0.02, "cpd00052": 0.02, "cpd00062": 0.02,
    "cpd00115": 0.01, "cpd00241": 0.01, "cpd00356": 0.01, "cpd00357": 0.01,
}


def parse_stoich(s):
    out = []
    for term in (s or "").split(";"):
        p = term.split(":")
        if len(p) < 3:
            continue
        try:
            c = float(p[0])
        except ValueError:
            continue
        if CPD_RE.fullmatch(p[1]):
            out.append((c, p[1], int(p[2]) if p[2].isdigit() else 0))
    return out


def build_universal(reactions_tsv, compounds_tsv, existing_ids):
    names = {}
    with open(compounds_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            names[row["id"]] = row.get("name", row["id"])
    uni = cobra.Model("universal")
    comp = {0: "c", 1: "e", 2: "p"}
    mets = {}
    with open(reactions_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("is_obsolete") or "0") == "1" or row["id"] in existing_ids:
                continue
            st = parse_stoich(row.get("stoichiometry"))
            if not st:
                continue
            r = Reaction(row["id"])
            rev = row.get("reversibility") or "="
            r.lower_bound = -1000.0 if rev in ("=", "<") else 0.0
            r.upper_bound = 1000.0 if rev in ("=", ">") else 0.0
            if r.lower_bound == 0 and r.upper_bound == 0:
                r.upper_bound = 1000.0
            sm = {}
            for c, cid, k in st:
                key = f"{cid}_{comp.get(k,'c')}"
                m = mets.get(key) or Metabolite(key, name=names.get(cid, cid), compartment=comp.get(k, "c"))
                mets[key] = m
                sm[m] = sm.get(m, 0) + c
            sm = {m: v for m, v in sm.items() if v != 0}
            if sm:
                r.add_metabolites(sm)
                uni.add_reactions([r])
    return uni


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT / "models/BOM_ss5_medium.xml"))
    ap.add_argument("--out", default=str(ROOT / "models/BOM_ss5_gapfilled.xml"))
    ap.add_argument("--reactions", default=str(ROOT / "refs/modelseed/reactions.tsv"))
    ap.add_argument("--compounds", default=str(ROOT / "refs/modelseed/compounds.tsv"))
    args = ap.parse_args()

    model = cobra.io.read_sbml_model(args.model)
    print(f"start: {len(model.reactions):,} reactions, {len(model.metabolites):,} metabolites")

    universal = build_universal(args.reactions, args.compounds, {r.id for r in model.reactions})
    print(f"universal pool: {len(universal.reactions):,} candidate reactions")

    added_total = set()
    for cid, label in TARGETS.items():
        key = f"{cid}_c"
        if key not in {m.id for m in model.metabolites}:
            model.add_metabolites([Metabolite(key, name=label, compartment="c")])
        met = model.metabolites.get_by_id(key)
        with model:
            dm = model.add_boundary(met, type="demand")
            model.objective = dm
            if model.slim_optimize() > 1e-6:
                print(f"  {label:<18} already producible")
                continue
        try:
            with model:
                dm = model.add_boundary(met, type="demand")
                model.objective = dm
                sol = gapfill(model, universal, demand_reactions=False, iterations=1)
            rxns = sol[0]
            print(f"  {label:<18} gapfilled with {len(rxns)} reaction(s): "
                  f"{', '.join(r.id for r in rxns[:6])}")
            for r in rxns:
                if r.id not in {x.id for x in model.reactions}:
                    c = r.copy()
                    c.notes["gapfilled"] = "true"
                    c.gene_reaction_rule = ""
                    model.add_reactions([c])
                    added_total.add(c.id)
        except Exception as e:
            print(f"  {label:<18} GAPFILL FAILED: {type(e).__name__}: {str(e)[:90]}")

    print(f"\ngapfilled reactions added: {len(added_total)}")

    # biomass objective from whatever is now reachable
    bio = Reaction("BIOMASS_fungal")
    bio.name = "Coarse fungal biomass (uncurated)"
    bio.lower_bound, bio.upper_bound = 0.0, 1000.0
    sm = {}
    missing = []
    for cid, coef in BIOMASS.items():
        key = f"{cid}_c"
        try:
            sm[model.metabolites.get_by_id(key)] = -coef
        except KeyError:
            missing.append(cid)
    bio.add_metabolites(sm)
    model.add_reactions([bio])
    model.objective = bio
    val = model.slim_optimize()
    print(f"biomass precursors in model: {len(sm)}/{len(BIOMASS)} (missing {len(missing)})")
    print(f"biomass flux on MNM medium: {val:.4f}" if val and val == val else "biomass flux: infeasible/zero")

    blocked = len(cobra.flux_analysis.find_blocked_reactions(model))
    print(f"blocked reactions: {blocked:,}/{len(model.reactions):,} "
          f"({100*blocked/len(model.reactions):.1f}%)")

    out = pathlib.Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(out))
    print(f"wrote {out}")
    print("NOTE: gapfilled reactions carry NO gene association and are tagged notes['gapfilled'].")


if __name__ == "__main__":
    sys.exit(main())
