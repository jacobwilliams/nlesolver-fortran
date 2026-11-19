# nlesolver-fortran technical note

`nlesolver_technote.tex` describes the theory behind the solver, its options, and their
performance, with recommended settings for common use cases.

## Contents

| Path | Description |
|---|---|
| `nlesolver_technote.tex` | the note (LaTeX) |
| `nlesolver_technote.pdf` | the compiled note |
| `benchmarks/nlesolver_benchmarks.f90` | benchmark program that produces all the data in the note |
| `data/` | benchmark results (CSV), read directly by the figures (pgfplots) |
| `scripts/make_tables.py` | generates the tables in `tables/` from the CSV files |
| `scripts/parabola.py` | generates the data for the parabola figure (`data/parabola_*.csv`) |
| `tables/` | generated LaTeX tables |

## Reproducing the results

From the repository root:

```sh
# 1. run the benchmarks (a few minutes; writes docs/data/*.csv)
fpm run --example nlesolver_benchmarks --profile release -- docs/data

# 2. regenerate the tables, and the parabola figure data
python3 docs/scripts/make_tables.py
python3 docs/scripts/parabola.py

# 3. build the PDF (needs a TeX distribution with pgfplots)
cd docs && latexmk -pdf nlesolver_technote.tex
```

Timings depend on the machine; the discussion in the note refers to the results in `data/`.
