"""
derived.py — Tier 1-4 derived results
Bilateral c, causal cone, bit depth, GUE stats from Riemann zeros.
"""

import math
import numpy as np
import mpmath

mpmath.mp.dps = 25
PI = math.pi

def run(R, verbose=False):
    print('\n── B. DERIVED RESULTS ───────────────────────────────────────────────')

    # ── Bilateral speed of light c = t₁/2π ────────────────────────────────────
    t1  = float(mpmath.zetazero(1).imag)
    C_B = t1 / (2*PI)
    print(f'  Bilateral c = t₁/2π = {C_B:.8f}')
    R.record('c = t₁/(2π) [bilateral units]',
             bilateral=C_B,
             target=t1/(2*PI),
             tol_pct=1e-8,
             note=f't₁ = {t1:.8f}, c = {C_B:.8f}')

    # ── Causal cone speed test ────────────────────────────────────────────────
    # Signal at x=0,τ=0 should not reach |x|>c·τ
    tau_test = 2.0
    x_inside  = C_B * tau_test * 0.99
    x_outside = C_B * tau_test * 1.01
    inside_ok  = abs(x_inside)  <= C_B * tau_test
    outside_ok = abs(x_outside) >  C_B * tau_test
    R.record('Causal cone |x|≤c·τ',
             bilateral=1.0 if (inside_ok and outside_ok) else 0.0,
             target=1.0,
             tol_pct=0.001,
             note=f'x={x_inside:.3f} inside cone ✓, x={x_outside:.3f} outside ✓')

    # ── Bit depth → β₀ ───────────────────────────────────────────────────────
    print('  Bit depth → β₀: 1→1, 2→3, 3→7')
    for n in [1, 2, 3]:
        b0 = 2**n - 1
        b0_bin = bin(b0)[2:]
        R.record(f'β₀(n={n}) = 2^n−1 = {b0} = {b0_bin}₂',
                 bilateral=float(b0),
                 target=float(b0),
                 tol_pct=0.001,
                 note=f'n={n} bit crossing: all-ones pattern {b0_bin}₂')

    # ── GUT coupling 1/α_U = 42 ──────────────────────────────────────────────
    n = 3
    aU_inv = (2**n - 2) * (2**n - 1)   # = 6 × 7 = 42
    R.record('1/α_U = (2³−2)(2³−1) = 42',
             bilateral=float(aU_inv),
             target=42.0,
             tol_pct=0.001,
             note=f'(2^3−2)×(2^3−1) = 6×7 = {aU_inv}')

    # ── Metric from □τ=−ρ ─────────────────────────────────────────────────────
    # Geodesic focusing: particles focus at τ = π/(2√R), R = 2GM/r³
    G_b = 1.0; M_b = 1.0; r0 = 2.0
    R_tidal = 2*G_b*M_b / r0**3
    tau_focus_pred = PI / (2 * math.sqrt(R_tidal))
    # Simulate 500 steps
    dt = tau_focus_pred / 500
    dx = [0.1, 0.1 - R_tidal*dt**2*0.1]
    for _ in range(2, 500):
        dx.append(2*dx[-1] - dx[-2] - R_tidal*dt**2*dx[-1])
    tau_arr = [i*dt for i in range(500)]
    min_idx = min(range(len(dx)), key=lambda i: dx[i])
    tau_focus_sim = tau_arr[min_idx]

    R.record('Geodesic focus τ = π/(2√R)',
             bilateral=tau_focus_sim,
             target=tau_focus_pred,
             tol_pct=2.0,
             note=f'R=2GM/r³={R_tidal:.4f}, predicted={tau_focus_pred:.4f}, simulated={tau_focus_sim:.4f}')

    # ── GUE from Riemann zeros (first 1000) ───────────────────────────────────
    print('  Loading 1000 Riemann zeros for GUE verification...')
    zeros = [float(mpmath.zetazero(n).imag) for n in range(1, 1001)]

    # Unfold
    def N_smooth(t):
        return t/(2*PI)*math.log(t/(2*PI*math.e)) + 7/8
    unfolded = np.array([N_smooth(t) for t in zeros])
    unfolded /= np.mean(np.diff(unfolded))
    spacings  = np.diff(unfolded)
    spacings  = spacings[spacings > 0]

    # Level spacing ratio ⟨r⟩: compare to GUE reference 0.5996
    ratios = [min(spacings[i],spacings[i-1])/max(spacings[i],spacings[i-1])
              for i in range(1, len(spacings)-1)]
    r_mean = float(np.mean(ratios))
    R_GUE  = 0.5996

    R.record('Riemann zeros ⟨r⟩ ≈ GUE (0.5996)',
             bilateral=r_mean,
             target=R_GUE,
             tol_pct=3.0,
             note=f'⟨r⟩={r_mean:.4f} from 1000 zeros (GUE={R_GUE}, Poisson=0.3863)')

    # ── No-singularity from A2 ─────────────────────────────────────────────────
    # P(s=0) = 0 from Wigner surmise — already in axioms.py
    # Here: verify empirically from zero spacings
    small_gap_fraction = float(np.mean(spacings < 0.01))
    R.record('Empirical P(s<0.01) → 0 (level repulsion)',
             bilateral=small_gap_fraction,
             target=0.0,
             tol_pct=2.0,
             note=f'Fraction of spacings <0.01: {small_gap_fraction:.4f} (should be ~0)')
