#!/usr/bin/env python3
"""
bilateral_verify.py
===================
Bilateral Mesh Framework — Full Verification Suite
Dunstan Low · ontologia.co.uk

Runs all derivations from three axioms to observables.
No free parameters. No fitting. Reports pass/fail with pulls.

Usage:
    python bilateral_verify.py
    python bilateral_verify.py --verbose
    python bilateral_verify.py --section rge
    python bilateral_verify.py --section all

Requirements:
    pip install mpmath scipy numpy sympy
"""

import sys
import math
import argparse
import importlib
import importlib.util
from pathlib import Path

# ── Colour output ─────────────────────────────────────────────────────────────
class C:
    PASS   = '\033[92m'   # green
    FAIL   = '\033[91m'   # red
    WARN   = '\033[93m'   # amber
    BOLD   = '\033[1m'
    BLUE   = '\033[94m'
    RESET  = '\033[0m'
    MUTED  = '\033[90m'

def col(text, code): return f"{code}{text}{C.RESET}"

# ── Import sub-modules ────────────────────────────────────────────────────────
def import_module(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / f"{name}.py")
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# ── Result collector ──────────────────────────────────────────────────────────
class Results:
    def __init__(self):
        self.passed  = []
        self.warned  = []
        self.failed  = []
        self.pending = []
        self.rows    = []   # (name, bilateral, target, dev_pct, pull, status, note)

    def record(self, name, bilateral, target, uncertainty=None,
               tol_pct=5.0, note='', pending=False):
        if pending:
            self.pending.append(name)
            self.rows.append((name, bilateral, target, None, None, 'PENDING', note))
            return
        if target is None:
            self.rows.append((name, bilateral, None, None, None, 'INFO', note))
            return
        dev = (bilateral - target) / (abs(target) if abs(target) > 1e-30 else 1.0) * 100
        pull = (bilateral - target) / uncertainty if uncertainty else None
        abs_pull = abs(pull) if pull else 0
        abs_dev  = abs(dev)

        if pull is not None:
            status = 'PASS' if abs_pull < 2 else ('WARN' if abs_pull < 3 else 'FAIL')
        else:
            status = 'PASS' if abs_dev < tol_pct else ('WARN' if abs_dev < 10 else 'FAIL')

        if status == 'PASS': self.passed.append(name)
        elif status == 'WARN': self.warned.append(name)
        else: self.failed.append(name)

        self.rows.append((name, bilateral, target, dev, pull, status, note))

    def print_table(self, verbose=False):
        print()
        w = [34, 14, 14, 10, 8, 8]
        hdr = ['Observable', 'Bilateral', 'Target', 'Dev %', 'Pull σ', 'Status']
        sep = '  '.join('─'*wi for wi in w)
        print(col('  '.join(h.ljust(w[i]) for i,h in enumerate(hdr)), C.BOLD))
        print(sep)
        for name, bil, tgt, dev, pull, status, note in self.rows:
            bil_s  = f'{bil:.8g}' if bil is not None else '—'
            tgt_s  = f'{tgt:.8g}' if tgt is not None else '—'
            dev_s  = f'{dev:+.4f}' if dev is not None else '—'
            pull_s = f'{pull:+.2f}' if pull is not None else '—'
            scol = C.PASS if status=='PASS' else C.FAIL if status=='FAIL' else C.WARN if status=='WARN' else C.MUTED
            row = '  '.join([
                name[:w[0]].ljust(w[0]),
                bil_s[:w[1]].ljust(w[1]),
                tgt_s[:w[2]].ljust(w[2]),
                dev_s[:w[3]].ljust(w[3]),
                pull_s[:w[4]].ljust(w[4]),
                col(status.ljust(w[5]), scol),
            ])
            print(row)
            if verbose and note:
                print(col(f'    {note}', C.MUTED))
        print(sep)

    def summary(self):
        total = len(self.passed)+len(self.warned)+len(self.failed)+len(self.pending)
        print()
        print(col('━'*70, C.BOLD))
        print(col('BILATERAL VERIFICATION SUMMARY', C.BOLD))
        print(col('━'*70, C.BOLD))
        print(f"  {col(str(len(self.passed)),  C.PASS)}  passed   (|pull| < 2σ or |dev| < 5%)")
        print(f"  {col(str(len(self.warned)),  C.WARN)}  warned   (2σ ≤ |pull| < 3σ)")
        print(f"  {col(str(len(self.failed)),  C.FAIL)}  failed   (|pull| ≥ 3σ)")
        print(f"  {col(str(len(self.pending)), C.MUTED)}  pending  (experiment not yet performed)")
        print(f"  {col(str(total),             C.BOLD)}  total")
        print()
        if not self.failed:
            print(col('  ✓ No contradictions. No free parameters. No fitting.', C.PASS))
        else:
            print(col(f'  ✗ {len(self.failed)} result(s) failed — see above', C.FAIL))
        print()

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Bilateral framework verification')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--section', '-s', default='all',
        choices=['all','axioms','derived','rge','observables','spectral',
                 'solitons','constraints','predictions'])
    args = parser.parse_args()

    print(col('━'*70, C.BOLD))
    print(col('  BILATERAL MESH FRAMEWORK — VERIFICATION SUITE', C.BOLD))
    print(col('  Dunstan Low · ontologia.co.uk', C.MUTED))
    print(col('━'*70, C.BOLD))
    print(col('  Three axioms. No free parameters. No fitting.\n', C.MUTED))

    R = Results()
    sections = [args.section] if args.section != 'all' else [
        'axioms','derived','rge','observables','spectral','solitons','constraints'
    ]

    for sec in sections:
        mod = import_module(sec)
        mod.run(R, verbose=args.verbose)

    R.print_table(verbose=args.verbose)
    R.summary()

if __name__ == '__main__':
    main()
