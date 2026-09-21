# Structural and Practical Identifiability of a Shared Metabolic Cancer ODE under Multi-Channel Noisy Observation Maps

**Thesis #9** (series label NP-01). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Can a shared metabolic ODE recover unique (or practically unique) parameters from multi-channel metabolomics-style outputs, or does identifiability collapse under realistic noisy maps?

The calculations use one shared four-state right-hand side. Schedule S records only the steady lactate/glucose ratio. Schedule M records glucose, pyruvate, lactate and a glutamine-linked pool at steady state, with correlated noise and a lineage loading. A short time course is included only as an upper bound. Observation noise is a synthetic surrogate, not a download of CCLE or DepMap.

On this surrogate the ratio has Fisher rank 1 of 5. The snapshot, after the loadings are profiled out, has rank 4 of 5. The missing direction is a common rescaling of the five rates, which the loading absorbs, so each coordinate profile is flat. Fixing one rate as a scale gauge is a different question, and so is a time course. Neither is a cell-line measurement.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md). An earlier sketch lives in [project-confluence](https://github.com/cloudynirvana/project-confluence/blob/main/docs/manuscript/structural_identifiability_ccle_manuscript.md). That draft is not the source of the ranks in Chapter Four.

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/identifiability.py` | Seeded Fisher and profile sketches (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Steady states, spectra, and profiles |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/identifiability.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the sketches. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Structural and practical identifiability of a shared metabolic cancer ODE under multi-channel noisy observation maps [Internet]. Thesis #9 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-09-ccle-metabolic-ode-identifiability

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
