# Bilateral Mesh Framework

**Three axioms. No free parameters. No fitting.**

A framework deriving Standard Model observables from three foundational axioms.
All results are computationally reproducible. The primary experimental test is JUNO (2031).

---

## Quick start

```bash
pip install -r requirements.txt
python bilateral_minimal.py    # 2 minutes — derives sin²θ_W from scratch
python bilateral_verify.py     # full verification suite
```

**Output of `bilateral_minimal.py`:**

```
============================================================
  BILATERAL MESH — MINIMAL VERIFICATION
  Deriving sin²θ_W from three axioms
============================================================

AXIOMS
  A1: Existence is relational  →  no preferred frame
  A2: No crossing is preferred →  U× = iσ_x (unique)
  A3: τ monotonically increasing → metric (−,+,+,+)

  ...

  Bilateral prediction:  sin²θ_W = 0.23122
  Observed (PDG 2024):   sin²θ_W = 0.23121 ± 4e-05
  Deviation:             +0.0043%
  Pull:                  +0.25σ

  Status: PASS ✓

  No free parameters. No fitting. Three axioms in.
```

---

## Three axioms

| Axiom | Statement | Consequence |
|-------|-----------|-------------|
| A1 | Existence is relational — no absolute position | Lorentz invariance; ⟨A⟩_γ = 2 exactly |
| A2 | No crossing is preferred — all crossings equivalent | GUE statistics; P(s=0)=0; no singularities |
| A3 | τ monotonically increasing — becoming-time forward only | Metric signature (−,+,+,+); causal cone |

---

## Confirmed predictions (PDG 2024)

| Observable | Bilateral | Observed | Pull |
|-----------|-----------|----------|------|
| sin²θ_W | 0.23122 | 0.23121 ± 0.00004 | +0.25σ |
| m_H (GeV) | 125.249 | 125.25 ± 0.17 | −0.01σ |
| K_l (Koide) | 2/3 | 0.666661 | 0.0009% |
| δ_CKM (rad) | arctan(5/2) = 1.1903 | 1.208 ± 0.058 | −0.31σ |
| \|V_us\| | 0.22537 | 0.22498 ± 0.00069 | +0.57σ |
| \|V_cb\| | 0.04221 | 0.04182 ± 0.00082 | +0.48σ |
| \|V_ub\| | 0.003724 | 0.003684 ± 0.00011 | +0.36σ |
| α_s(M_Z) | 0.11658 | 0.1179 ± 0.0010 | −1.32σ ¹ |

¹ One-loop only. Two-loop β₁ is an open calculation.

**Zero free parameters. Zero fitting. Zero contradictions.**

---

## Pending predictions

| Observable | Bilateral prediction | Experiment | Timeline |
|-----------|---------------------|------------|----------|
| Neutrino ordering | Inverted ordering, m₃ = 0 exactly | JUNO / Hyper-K | ~2027–2031 |
| m₃ (lightest ν mass) | 0.000 eV exactly, m₁ ≈ m₂ ≈ 49.5 meV | Planck / CMB-S4 | ~2030s |

**Falsification:** Normal neutrino ordering confirmed at > 3σ falsifies the framework.

---

## Key derivations

```
β₀(SU3) = 2³ − 1         = 7 = 111₂   ← all-ones 3-bit pattern (A2)
β₀(SU2) = 2² − 1         = 3 = 11₂    ← all-ones 2-bit pattern (A2)
1/α_U   = (2³−2)(2³−1)   = 42          ← bilateral crossing product (A2)
1/α_s   = 42 − 7×30/(2π) = 8.577       ← prime rung k=30, p=113
1/α_2   = 42 − 3×25/(2π) = 30.06       ← prime rung k=25, p=97
sin²θ_W = α₂/(α_s+α₂)   = 0.23122     ← +0.25σ from PDG
m_H     = 125.000 + 0.499 = 125.249 GeV ← 0.001% from LHC
c       = t₁/(2π)         = 2.2496...  ← first Riemann zero / 2π
```

---

## Files

| File | Purpose |
|------|---------|
| `bilateral_minimal.py` | **Start here.** Derives sin²θ_W in ~80 lines. |
| `bilateral_spec.py` | Formal spec for independent reimplementation. Self-verifying. |
| `bilateral_verify.py` | Full verification suite runner. |
| `axioms.py` | A1 (Lorentz), A2 (GUE, P(s=0)=0), A3 (metric signature) |
| `derived.py` | Crossing operator, causal cone, bit depth, geodesic focusing |
| `rge.py` | β₀, 1/α_U, RGE, prime ladder, asymptotic freedom |
| `observables.py` | All 8 predictions vs PDG 2024 with pulls |
| `spectral.py` | Riemann zero statistics, rigidity, selection rule |
| `solitons.py` | Kink rest energy, topological charge, energy conservation |
| `constraints.py` | Parameter space scan (0/3120 continuous combinations pass) |

---

## Independent reimplementation

`bilateral_spec.py` is a language-neutral specification. It states the three axioms,
the five derivation steps, and the eight verification targets as a readable document.
Reimplement in any language (Julia, C++, Mathematica) and verify against the targets.

```bash
python bilateral_spec.py    # verifies the spec is self-consistent (13/13 checks)
```

---

## Run options

```bash
# Quick demo — derives sin²θ_W and Higgs mass from axioms
python bilateral_minimal.py

# Full suite
python bilateral_verify.py

# Single section
python bilateral_verify.py --section observables
python bilateral_verify.py --section rge
python bilateral_verify.py --section axioms

# With derivation notes
python bilateral_verify.py --verbose

# Available sections: axioms, derived, rge, observables, spectral, solitons, constraints
```

---

## Website

Simulations, visualisations, and papers: [ontologia.co.uk](https://ontologia.co.uk)

Includes interactive pages for vacuum fluctuations, solitons, spacetime, RG flow,
horizons, constraint search, SM interactions, cosmological timeline, proof trace,
perturbation engine, prediction interface, and the zeta 720° diagram.

---

## Author

Dunstan Low · independent scholar · [ontologia.co.uk](https://ontologia.co.uk)

---

## Licence

MIT
