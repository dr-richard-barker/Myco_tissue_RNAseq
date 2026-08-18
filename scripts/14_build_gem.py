#!/usr/bin/env python3
"""Phase 6 -- draft genome-scale metabolic reconstruction for Pleurotus ostreatus PC9.15.

No published P. ostreatus GEM exists (literature search found none), so this is a de novo
draft. The intended route was ModelSEEDpy/eggNOG, but eggnogdb.embl.de is unreachable from
this machine, so the reconstruction is assembled directly:

    proteome --DIAMOND--> Swiss-Prot --EC--> ModelSEED reactions --> COBRA model

Gene-protein-reaction rules are built from the EC mapping, so every reaction is traceable to
the PC9.15 proteins that justify it. This is a DRAFT: it is unfitted, ungapfilled and carries
no organism-specific biomass objective, all of which need curation before flux predictions
are trustworthy. It is the scaffold that Phase 6 contextualisation will sit on, not a
finished model.
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

import cobra

ROOT = pathlib.Path(__file__).resolve().parents[1]
CPD_RE = re.compile(r"cpd\d+")


def load_ec_to_genes(func_tsv):
    ec2genes = collections.defaultdict(set)
    with open(func_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row["EC"]:
                continue
            for ec in row["EC"].split(";"):
                ec = ec.strip()
                # skip incomplete EC classes such as 1.1.1.- : they map to hundreds of
                # unrelated reactions and would flood the draft with false positives
                if ec and "-" not in ec:
                    ec2genes[ec].add(row["protein_id"])
    return ec2genes


def load_compounds(path):
    names, formulas, charges = {}, {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = row["id"]
            names[cid] = row.get("name", cid)
            f = row.get("formula", "")
            formulas[cid] = "" if f in ("null", "None") else f
            try:
                charges[cid] = int(float(row.get("charge") or 0))
            except ValueError:
                charges[cid] = 0
    return names, formulas, charges


def parse_stoich(stoich):
    """ModelSEED 'n:cpdXXXXX:compartment:...' entries, semicolon separated."""
    out = []
    for term in stoich.split(";"):
        parts = term.split(":")
        if len(parts) < 3:
            continue
        try:
            coef = float(parts[0])
        except ValueError:
            continue
        cid, comp = parts[1], parts[2]
        if not CPD_RE.fullmatch(cid):
            continue
        out.append((coef, cid, int(comp) if comp.isdigit() else 0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--functional", default=str(ROOT / "results/annotation/pc915_functional.tsv"))
    ap.add_argument("--reactions", default=str(ROOT / "refs/modelseed/reactions.tsv"))
    ap.add_argument("--compounds", default=str(ROOT / "refs/modelseed/compounds.tsv"))
    ap.add_argument("--out", default=str(ROOT / "models/PC9.15_draft.xml"))
    args = ap.parse_args()

    ec2genes = load_ec_to_genes(args.functional)
    print(f"complete EC numbers with gene support: {len(ec2genes):,}")

    names, formulas, charges = load_compounds(args.compounds)
    print(f"ModelSEED compounds loaded: {len(names):,}")

    model = cobra.Model("PC9_15_draft")
    comp_names = {0: "c", 1: "e", 2: "p"}
    metabolites = {}
    added, skipped_obsolete, skipped_nostoich = 0, 0, 0
    rxn_genes = {}

    with open(args.reactions) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ecs = (row.get("ec_numbers") or "").strip()
            if not ecs or ecs in ("null", "None"):
                continue
            hit = {e.strip() for e in ecs.split("|")} | {e.strip() for e in ecs.split(";")}
            genes = set()
            matched = []
            for e in hit:
                if e in ec2genes:
                    genes |= ec2genes[e]
                    matched.append(e)
            if not genes:
                continue
            if (row.get("is_obsolete") or "0") == "1":
                skipped_obsolete += 1
                continue
            st = parse_stoich(row.get("stoichiometry") or "")
            if not st:
                skipped_nostoich += 1
                continue

            rid = row["id"]
            rxn = cobra.Reaction(rid)
            rxn.name = (row.get("name") or rid)[:200]
            rev = (row.get("reversibility") or "=")
            rxn.lower_bound = -1000.0 if rev in ("=", "<") else 0.0
            rxn.upper_bound = 1000.0 if rev in ("=", ">") else 0.0
            if rxn.lower_bound == 0.0 and rxn.upper_bound == 0.0:
                rxn.upper_bound = 1000.0

            stoich_map = {}
            for coef, cid, comp in st:
                c = comp_names.get(comp, "c")
                key = f"{cid}_{c}"
                met = metabolites.get(key)
                if met is None:
                    met = cobra.Metabolite(key, name=names.get(cid, cid), compartment=c)
                    if formulas.get(cid):
                        met.formula = formulas[cid]
                    met.charge = charges.get(cid, 0)
                    metabolites[key] = met
                stoich_map[met] = stoich_map.get(met, 0) + coef
            stoich_map = {m: v for m, v in stoich_map.items() if v != 0}
            if not stoich_map:
                skipped_nostoich += 1
                continue
            rxn.add_metabolites(stoich_map)
            model.add_reactions([rxn])
            rxn_genes[rid] = genes
            added += 1

    # GPR rules: any of the proteins carrying the matching EC can support the reaction
    for rid, genes in rxn_genes.items():
        model.reactions.get_by_id(rid).gene_reaction_rule = " or ".join(sorted(genes))

    # exchanges for everything in the extracellular compartment, so the draft is at least
    # simulatable; the medium is not yet constrained to the real culture substrate
    n_ex = 0
    for met in list(model.metabolites):
        if met.compartment == "e":
            model.add_boundary(met, type="exchange")
            n_ex += 1

    print(f"reactions added        : {added:,}")
    print(f"  skipped (obsolete)   : {skipped_obsolete:,}")
    print(f"  skipped (no stoich)  : {skipped_nostoich:,}")
    print(f"metabolites            : {len(model.metabolites):,}")
    print(f"genes                  : {len(model.genes):,}")
    print(f"exchange reactions     : {n_ex:,}")

    compartments = collections.Counter(m.compartment for m in model.metabolites)
    print(f"compartments           : {dict(compartments)}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(out))
    print(f"wrote {out}")

    # connectivity sanity check: a draft with no blocked-free core is not worth carrying on
    orphan = [m.id for m in model.metabolites if len(m.reactions) == 1]
    print(f"metabolites in exactly one reaction (dead ends): {len(orphan):,} "
          f"({100 * len(orphan) / max(len(model.metabolites), 1):.1f}%)")


if __name__ == "__main__":
    sys.exit(main())
