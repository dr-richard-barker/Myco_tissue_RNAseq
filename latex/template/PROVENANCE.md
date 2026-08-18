# Springer Nature LaTeX template — provenance

| | |
|---|---|
| Source | https://www.springernature.com/gp/authors/campaigns/latex-author-support |
| Archive | `springernature.zip` (December 2024 version, 901,814 bytes) |
| Retrieved | 2026-08-18 |
| Class | `sn-jnl.cls` v0.1 |
| Licence | LaTeX Project Public License (LPPL) 1.3c or later |

## Class options used

```latex
\documentclass[sn-nature,pdflatex]{sn-jnl}
```

- **`sn-nature`** — the Nature Portfolio reference style, which is the correct choice for
  npj titles. Bibliography style `bst/sn-nature.bst`.
- **`pdflatex`** — **required here, and not for the reason the name suggests.** Without it the
  class executes `\RequirePackage[hyphenbreaks]{breakurl}`, a dvips-only package that fails
  under any modern engine with `Undefined control sequence \headerps@out`. The option gates
  that `\RequirePackage` off. It means "a PDF-producing engine", which XeTeX satisfies.
- **`lineno`** — added to review builds only.

## Build engine

`tectonic` 0.17.0 (conda-forge, `env-tex`). Self-contained; fetches LaTeX packages on demand,
avoiding a full TeX Live installation.

## Two constraints that follow from using XeTeX

1. **Figures must be PDF, not EPS.** The shipped `sn-article.tex` example fails at
   `image inclusion failed for "fig.eps"` because xdvipdfmx cannot embed EPS without
   Ghostscript. All figures here are generated as vector PDF.
2. **The preamble must load `amsmath`** (and `graphicx`, `booktabs`). The class does not
   load them itself, and a minimal document fails with an undefined `\allowdisplaybreaks`.
