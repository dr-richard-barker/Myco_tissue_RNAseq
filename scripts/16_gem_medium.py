#!/usr/bin/env python3
"""Phase 6 -- put the draft GEM on the real culture medium and measure what it can make.

The medium is the mycoponic nutrient medium (MNM v3) of Porterfield et al. 2026
(doi:10.1002/biot.70184), used to grow these cultures -- see metadata/culture_conditions.md.
It is a complex rich medium (corn syrup, malt extract, peptone, tryptic soy broth, ammonium
sulfate, gypsum, enzymatically pre-digested cellulose and oak sawdust), so it is modelled as
permitting uptake of common sugars, the 20 proteinogenic amino acids, ammonium, sulfate,
phosphate and mineral ions rather than a single defined carbon source.

The draft from 14_build_gem.py is 81.5% blocked, overwhelmingly because ModelSEED transport
reactions carry no EC number and so were never picked up by the EC-driven build. This script
adds transport and exchange reactions for the medium components and for metabolites the
network already contains, then reports:

  * how much of the network becomes able to carry flux
  * which standard biomass precursors the model can synthesise from the medium

It does NOT claim a validated growth prediction. There is still no curated biomass objective
and no mass-balance curation; the precursor test is a connectivity diagnostic, and is
reported as such.
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

import cobra
from cobra import Metabolite, Reaction

ROOT = pathlib.Path(__file__).resolve().parents[1]
CPD_RE = re.compile(r"cpd\d+")

# ModelSEED compound ids for the medium and for the precursors we test.
MEDIUM = {
    "cpd00027": "D-glucose", "cpd00082": "D-fructose", "cpd00179": "maltose",
    "cpd00076": "sucrose", "cpd00013": "ammonium", "cpd00048": "sulfate",
    "cpd00009": "phosphate", "cpd00001": "water", "cpd00007": "O2",
    "cpd00011": "CO2", "cpd00067": "H+", "cpd00099": "chloride",
    "cpd00063": "Ca2+", "cpd00205": "K+", "cpd00254": "Mg2+", "cpd00971": "Na+",
    "cpd10515": "Fe2+", "cpd00030": "Mn2+", "cpd00034": "Zn2+", "cpd00058": "Cu2+",
    "cpd00149": "Co2+", "cpd00220": "riboflavin", "cpd00305": "thiamine",
    # Biotin is supplied, not synthesised. It stays unreachable after targeted gapfilling,
    # and that is biologically right rather than a modelling gap: many fungi are biotin
    # auxotrophs, and MNM v3 contains malt extract, peptone and tryptic soy broth, all of
    # which supply B-vitamins. Adding a biosynthesis route would have been the wrong fix.
    "cpd00104": "biotin",
    # peptone / tryptic soy broth supply free amino acids
    "cpd00035": "L-alanine", "cpd00051": "L-arginine", "cpd00132": "L-asparagine",
    "cpd00041": "L-aspartate", "cpd00084": "L-cysteine", "cpd00053": "L-glutamine",
    "cpd00023": "L-glutamate", "cpd00033": "glycine", "cpd00119": "L-histidine",
    "cpd00322": "L-isoleucine", "cpd00107": "L-leucine", "cpd00039": "L-lysine",
    "cpd00060": "L-methionine", "cpd00066": "L-phenylalanine", "cpd00129": "L-proline",
    "cpd00054": "L-serine", "cpd00161": "L-threonine", "cpd00065": "L-tryptophan",
    "cpd00069": "L-tyrosine", "cpd00156": "L-valine",
}

# Standard biomass precursors -- can the network reach them from the medium?
PRECURSORS = {
    "cpd00002": "ATP", "cpd00038": "GTP", "cpd00052": "CTP", "cpd00062": "UTP",
    "cpd00115": "dATP", "cpd00241": "dGTP", "cpd00356": "dCTP", "cpd00357": "dTTP",
    "cpd00035": "L-alanine", "cpd00023": "L-glutamate", "cpd00033": "glycine",
    "cpd00054": "L-serine", "cpd00161": "L-threonine", "cpd00129": "L-proline",
    "cpd00066": "L-phenylalanine", "cpd00069": "L-tyrosine", "cpd00065": "L-tryptophan",
    "cpd00060": "L-methionine", "cpd00084": "L-cysteine", "cpd00039": "L-lysine",
    "cpd00051": "L-arginine", "cpd00119": "L-histidine", "cpd00322": "L-isoleucine",
    "cpd00107": "L-leucine", "cpd00156": "L-valine", "cpd00041": "L-aspartate",
    "cpd00053": "L-glutamine", "cpd00132": "L-asparagine",
    "cpd00010": "CoA", "cpd00003": "NAD", "cpd00006": "NADP", "cpd00015": "FAD",
    "cpd00037": "UDP-GlcNAc", "cpd00104": "biotin",
}


def parse_stoich(stoich):
    out = []
    for term in (stoich or "").split(";"):
        parts = term.split(":")
        if len(parts) < 3:
            continue
        try:
            coef = float(parts[0])
        except ValueError:
            continue
        cid, comp = parts[1], parts[2]
        if CPD_RE.fullmatch(cid):
            out.append((coef, cid, int(comp) if comp.isdigit() else 0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT / "models/PC9.15_draft.xml"))
    ap.add_argument("--reactions", default=str(ROOT / "refs/modelseed/reactions.tsv"))
    ap.add_argument("--compounds", default=str(ROOT / "refs/modelseed/compounds.tsv"))
    ap.add_argument("--out", default=str(ROOT / "models/PC9.15_medium.xml"))
    args = ap.parse_args()

    model = cobra.io.read_sbml_model(args.model)
    before_rxn = len(model.reactions)
    blocked_before = len(cobra.flux_analysis.find_blocked_reactions(model))
    print(f"draft: {before_rxn:,} reactions, {blocked_before:,} blocked "
          f"({100 * blocked_before / before_rxn:.1f}%)")

    names = {}
    with open(args.compounds) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            names[row["id"]] = row.get("name", row["id"])

    comp_names = {0: "c", 1: "e", 2: "p"}
    present = {m.id for m in model.metabolites}

    def get_met(cid, c):
        key = f"{cid}_{c}"
        if key in present:
            return model.metabolites.get_by_id(key)
        met = Metabolite(key, name=names.get(cid, cid), compartment=c)
        model.add_metabolites([met])
        present.add(key)
        return met

    # ---- add ModelSEED transport reactions touching metabolites we already have ----
    added_t = 0
    with open(args.reactions) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("is_transport") or "0") != "1":
                continue
            if (row.get("is_obsolete") or "0") == "1":
                continue
            st = parse_stoich(row.get("stoichiometry"))
            if not st:
                continue
            cids = {cid for _c, cid, _k in st}
            # only transporters whose cytosolic partner the network already knows about,
            # so we extend connectivity rather than bolting on unrelated chemistry
            if not any(f"{cid}_c" in present for cid in cids):
                continue
            rid = row["id"]
            if rid in model.reactions:
                continue
            rxn = Reaction(rid)
            rxn.name = (row.get("name") or rid)[:200]
            rev = row.get("reversibility") or "="
            rxn.lower_bound = -1000.0 if rev in ("=", "<") else 0.0
            rxn.upper_bound = 1000.0 if rev in ("=", ">") else 0.0
            if rxn.lower_bound == 0.0 and rxn.upper_bound == 0.0:
                rxn.upper_bound = 1000.0
            sm = {}
            for coef, cid, comp in st:
                met = get_met(cid, comp_names.get(comp, "c"))
                sm[met] = sm.get(met, 0) + coef
            sm = {m: v for m, v in sm.items() if v != 0}
            if not sm:
                continue
            rxn.add_metabolites(sm)
            model.add_reactions([rxn])
            added_t += 1
    print(f"transport reactions added: {added_t:,}")

    # ---- medium: exchanges open for MNM components, closed elsewhere ----
    existing_ex = {r.id for r in model.reactions if r.boundary}
    opened = []
    for cid, label in MEDIUM.items():
        met = get_met(cid, "e")
        exid = f"EX_{met.id}"
        if exid not in existing_ex:
            try:
                model.add_boundary(met, type="exchange")
            except ValueError:
                continue
        model.reactions.get_by_id(exid).lower_bound = -10.0
        opened.append(label)

    for r in model.reactions:
        if r.boundary and r.id not in {f"EX_{c}_e" for c in MEDIUM}:
            r.lower_bound = 0.0   # no uptake of anything not in the medium
    print(f"medium exchanges opened: {len(opened)} (MNM v3, Porterfield et al. 2026)")

    blocked_after = len(cobra.flux_analysis.find_blocked_reactions(model))
    print(f"after medium+transport: {len(model.reactions):,} reactions, "
          f"{blocked_after:,} blocked ({100 * blocked_after / len(model.reactions):.1f}%)")
    print(f"reactions able to carry flux: {before_rxn - blocked_before:,} -> "
          f"{len(model.reactions) - blocked_after:,}")

    # ---- precursor reachability: add a temporary demand and see if it can carry flux ----
    reachable, unreachable = [], []
    for cid, label in PRECURSORS.items():
        key = f"{cid}_c"
        if key not in present:
            unreachable.append((label, "not in network"))
            continue
        met = model.metabolites.get_by_id(key)
        with model:
            dm = model.add_boundary(met, type="demand")
            model.objective = dm
            try:
                val = model.slim_optimize()
            except Exception:
                val = 0.0
        if val and val > 1e-6:
            reachable.append(label)
        else:
            unreachable.append((label, "blocked"))

    print(f"\nbiomass precursors synthesisable from MNM: "
          f"{len(reachable)}/{len(PRECURSORS)}")
    print("  reachable  :", ", ".join(sorted(reachable)) or "(none)")
    missing = collections.Counter(r for _l, r in unreachable)
    print(f"  unreachable: {len(unreachable)} ({dict(missing)})")
    for label, why in sorted(unreachable)[:15]:
        print(f"     {label:<22} {why}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(out))
    print(f"\nwrote {out}")
    print("NOTE: connectivity diagnostic only -- no curated biomass objective, no "
          "mass-balance curation. Not a validated growth prediction.")


if __name__ == "__main__":
    sys.exit(main())
