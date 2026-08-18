# Culture conditions

Source: Porterfield et al. (2026), *Mycoponics: Controlled Bioproduction Utilizing
Biophysical, Solid-State, Liquid Nutrient Delivery*, Biotechnology Journal,
[10.1002/biot.70184](https://doi.org/10.1002/biot.70184) — open access at
[PMC12878557](https://pmc.ncbi.nlm.nih.gov/articles/PMC12878557/).

## Strain

***Pleurotus ostreatus* cv. "Harbor Blue P01"** — a commercial **blue oyster** cultivar.
(The paper also used "Snow-White Oyster Po3"; this dataset is the blue oyster.)

This matters for the reference choice. `BOM_ss5` / `BOM_ss14` (`GCA_056149245.1` /
`GCA_056149315.1`) are Academia Sinica **B**lue **O**yster **M**ushroom single-spore
isolates, and they mapped 3.5 points better than PC9.15 on raw unique rate while scoring
within 0.02 points of each other — the signature of the two nuclei of one dikaryon.
PC9.15 was selected only because its annotation carries real UTRs, which dominated the gene
assignment comparison. See NOTES.md §4 and the re-test with matched annotations.

## Mycoponic nutrient medium (MNM v3)

Per 1500 mL water:

| Component | Amount | Provides |
|---|---|---|
| Corn syrup | 110.0 g | glucose, maltose, fructose |
| Malt extract | 8.0 g | maltose, maltodextrin, amino acids |
| Peptone | 7.0 g | amino acids, peptides |
| Tryptic soy broth | 7.0 g | amino acids, glucose, salts |
| Microcrystalline cellulose | 5.0 g | glucose (enzymatically pre-digested) |
| Oak sawdust | 5.0 g | cellulose, hemicellulose/xylan, lignin |
| Nutritional enzymes | 2.0 g | hydrolytic consortium |
| Ammonium sulfate | 0.5 g | NH4+, SO4(2-) |
| Gypsum (CaSO4) | 0.7 g | Ca(2+), SO4(2-) |

Cellulose and oak sawdust were pre-incubated with the hydrolytic enzyme consortium (2.3 g)
at 45 °C for 24 h before mixing. All media autoclaved before introduction to the tubes.

**This is a complex/undefined rich medium**, not a minimal defined one. For gapfilling, the
GEM medium is therefore modelled as permitting uptake of: mono- and disaccharides (glucose,
fructose, maltose), the 20 proteinogenic amino acids, ammonium, sulfate, phosphate, and
standard mineral/trace ions — rather than a single carbon source. Lignocellulose components
(cellulose, xylan, lignin) are available but require secreted CAZymes and oxidative enzymes,
which is directly relevant to the tissue biology.

## Physical setup

- Ceramic MiniTubes, 10 cm length x 5 cm diameter, pore size < 300 nm
- Filled 50% by volume with granular activated carbon, sealed with 3D-printed caps
- Air-phase/solid-state interface; "persistent-filtration-defense" excludes microbes
- **16 °C**, **85% relative humidity**
- **CO2 <= 1000 ppm**, HEPA filtration
- Cultures viable over 7 months without contamination

## Implications for the analysis

- Growth at 16 °C and <=1000 ppm CO2 is cool and well-aerated; not hypoxic.
- The activated carbon adsorbs organics, which may shape what the mycelium actually
  encounters versus the nominal recipe.
- Oak sawdust plus microcrystalline cellulose means lignocellulolytic machinery (CAZymes,
  laccases, Mn/versatile peroxidases) is expected to be active and is a priority annotation
  target for the per-tissue models.
- The exudate tissues (exuding mycelium, exudophore) are sampled from a system with a liquid
  nutrient interface, so secretome and transport functions are the biology of interest.
