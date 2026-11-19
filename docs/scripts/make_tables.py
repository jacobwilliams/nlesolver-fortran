#!/usr/bin/env python3
"""Generate the LaTeX tables of the technical note from the benchmark CSV files.

Usage (from the repository root, after running the benchmarks):

    python3 docs/scripts/make_tables.py

Reads docs/data/*.csv and writes docs/tables/*.tex.
"""

import csv
import os
import platform
import subprocess

DOCS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
DATA = os.path.join(DOCS, 'data')
TABLES = os.path.join(DOCS, 'tables')

NAMES = {'dense': 'dense', 'lusol': 'LUSOL', 'lsqr': 'LSQR', 'lsmr': 'LSMR',
         'craig': 'CRAIG', 'cholesky': 'Cholesky'}
STATUS = {1: 'converged', 2: 'no progress', 3: 'max iter', 5: 'uphill'}


def read(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def sci(v, digits=1):
    """number in LaTeX scientific notation"""
    v = float(v)
    if v == 0.0:
        return '0'
    s = f'{v:.{digits}e}'
    mant, exp = s.split('e')
    return f'${mant}\\times10^{{{int(exp)}}}$'


def sec(v):
    v = float(v)
    if v < 0.01:
        return f'{v:.4f}'
    if v < 1:
        return f'{v:.3f}'
    if v < 100:
        return f'{v:.2f}'
    return f'{v:.0f}'


def status(s):
    s = int(s)
    return STATUS.get(s, f'error ({s})')


def write(name, text):
    os.makedirs(TABLES, exist_ok=True)
    with open(os.path.join(TABLES, name), 'w') as f:
        f.write(text)


def setup():
    try:
        cpu = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'],
                             capture_output=True, text=True).stdout.strip()
    except OSError:
        cpu = platform.processor()
    try:
        fc = subprocess.run(['gfortran', '--version'], capture_output=True,
                            text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        fc = 'gfortran'
    write('setup.tex',
          '\\paragraph{Setup.} All timings are wall-clock times on a single core of an '
          f'{cpu} ({platform.system()} {platform.machine()}), compiled with {fc.replace("_", chr(92)+"_")} '
          'using the fpm \\code{release} profile, in double precision. Each time is the best of '
          'three runs for $n \\le 50{,}000$ and of one run otherwise. The benchmark program is '
          '\\code{docs/benchmarks/nlesolver\\_benchmarks.f90}; see \\code{docs/README.md} to '
          'reproduce the results.\n')


def scaling_table():
    rows = [r for r in read('scaling_detail.csv') if int(r['n']) == 50000]
    out = ['\\begin{table}[htbp]', '\\centering',
           '\\caption{Banded underdetermined problem, $n = 50{,}000$, $m = 30{,}000$ '
           '(backtracking line search). $\\norm{x - x_0}$ is the distance from the starting point '
           'to the computed solution.}',
           '\\label{tab:scaling}', '\\small',
           '\\begin{tabular}{@{}llrrrrrr@{}}', '\\toprule',
           'solver & status & iterations & $F$ evals & $J$ evals & $\\norm{F}$ & $\\norm{x-x_0}$ & time [s] \\\\',
           '\\midrule']
    for r in rows:
        out.append(f"{NAMES[r['solver']]} & {status(r['istat'])} & {r['iters']} & {r['fevals']} & "
                   f"{r['jevals']} & {sci(r['fnorm'])} & {float(r['dist']):.4f} & {sec(r['time'])} \\\\")
    out += ['\\bottomrule', '\\end{tabular}', '\\end{table}', '']
    write('scaling_table.tex', '\n'.join(out))


def steps_table():
    rows = read('steps_detail.csv')
    sm = {'1': 'constant', '2': 'backtracking', '3': 'exact', '4': 'fixed points'}
    out = ['\\begin{table}[htbp]', '\\centering',
           '\\caption{Step-length strategies with and without Schubert updates, banded problem with '
           '$n = 50{,}000$.}',
           '\\label{tab:steps}', '\\small',
           '\\begin{tabular}{@{}lllrrrrr@{}}', '\\toprule',
           'solver & step mode & updates & status & iterations & $F$ evals & $J$ evals & time [s] \\\\',
           '\\midrule']
    last = None
    for r in rows:
        if last is not None and r['solver'] != last:
            out.append('\\midrule')
        last = r['solver']
        out.append(f"{NAMES[r['solver']]} & {sm[r['step_mode']]} & {'yes' if r['broyden']=='1' else 'no'} & "
                   f"{status(r['istat'])} & {r['iters']} & {r['fevals']} & {r['jevals']} & {sec(r['time'])} \\\\")
    out += ['\\bottomrule', '\\end{tabular}', '\\end{table}', '']
    write('steps_table.tex', '\n'.join(out))


def sparsity_table():
    rows = read('sparsity.csv')
    out = ['\\begin{table}[htbp]', '\\centering',
           '\\caption{Effect of the sparsity pattern on the Cholesky solver: bandwidth of $JJ^T$ in the '
           'given order and after RCM, envelope size of $L$, and total solve times (Cholesky and CRAIG).}',
           '\\label{tab:sparsity}', '\\small',
           '\\begin{tabular}{@{}lrrrrrrr@{}}', '\\toprule',
           'problem & $n$ & $m$ & bandwidth & after RCM & envelope & Cholesky [s] & CRAIG [s] \\\\',
           '\\midrule']
    for r in rows:
        name = r['case'] if r['case'] != 'grid' else f"grid, $k={r['param']}$"
        out.append(f"{name} & {int(r['n']):,} & {int(r['m']):,} & {int(r['bw_natural']):,} & "
                   f"{int(r['bw_rcm']):,} & {int(r['envelope']):,} & {sec(r['chol_time'])} & "
                   f"{sec(r['craig_time'])} \\\\".replace(',', '{,}'))
    out += ['\\bottomrule', '\\end{tabular}', '\\end{table}', '']
    write('sparsity_table.tex', '\n'.join(out))


def envelope_polygons():
    """outline of the envelope (lower triangle, after RCM) for the sparsity pattern figures"""
    for name in ('shuffled', 'grid'):
        rows = read(f'spy_{name}_env.csv')
        pts = []
        for r in rows:                        # left edge: first column of each row
            i, f = int(r['row']), int(r['first'])
            pts += [(f - 0.5, i - 0.5), (f - 0.5, i + 0.5)]
        for r in reversed(rows):              # right edge: the diagonal
            i = int(r['row'])
            pts += [(i + 0.5, i + 0.5), (i + 0.5, i - 0.5)]
        with open(os.path.join(DATA, f'spy_{name}_env_poly.csv'), 'w') as f:
            f.write('x,y\n')
            for x, y in pts:
                f.write(f'{x},{y}\n')


if __name__ == '__main__':
    envelope_polygons()
    setup()
    scaling_table()
    steps_table()
    sparsity_table()
    print('tables written to', os.path.normpath(TABLES))
