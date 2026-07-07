# CR Two-Moment vs Athena++: Parity Plan

Follow-up to the pre-PR review of the CR_STREAMING module. Each known GAMER-vs-Athena
difference is classified by **difficulty**, **impact on matching**, and a **recommendation**:
`FIX` (do it), `FIX-LATER` (only needed for live-gas/coupled parity), or
`ACCEPT` (structural GAMER-framework difference; ask for sign-off, like the B.C. cadence).

Item numbers (B1..B8, C1..C7) refer to the review report.

## Summary table

| # | Item | Difficulty | When it matters | Recommendation |
|---|------|-----------|-----------------|----------------|
| B1 | TINY_NUMBER: DBL_MIN vs 1e-20 | Easy | Always (incl. current static tests) | **FIX first** |
| C1 | `va<=TINY` fallback + extra `Ec>TINY` guard | Easy | Edge cells (floored Ec, B~0) | **FIX with B1** |
| B7a | ADV_* passive-advection fluxes (HLLE skip commented out) | Easy | Always (stored fields, patch edges) | **FIX second** |
| B2 | Gas state fed to source solve (t^n / half-step vs stage-updated) | Easy–Moderate | Live gas only (CR_SOURCE=1 or moving gas) | FIX-LATER (bundle 1) |
| B3 | Missing gas back-reaction at half step | Easy–Moderate | Live gas + CR_SOURCE=1 | FIX-LATER (bundle 1) |
| B4 | ec_source: central-diff vs flux-divergence grad(Pc) | Moderate | Multi-D, moving gas, oblique B | FIX-LATER (bundle 1) |
| B5 | Half-step source uses half-step B vs Athena t^n B | Moderate (buffer aliasing) | Evolving B only | FIX-LATER (bundle 2) |
| B6 | Floor vs gas-energy update ordering in source | Trivial | Only when source drives Ec<0 with CR_SOURCE=1 | Decide with boss (match = trivial; GAMER's version conserves energy better) |
| B8 | MinMod retry reuses mutated sigma; 1st-order-flux-corr fallback | Easy (guard) / Moderate (retry) | Solver-failure paths only | Guard now, document retry |
| B7b | Stale ADV_* in outermost ghost ring | Hard (structural) | Patch-boundary cells, every step | **ACCEPT** (same class as B.C. cadence) |
| C2 | dt criterion semantics | Moderate (framework) | Only if gas speed ~ vmax | ACCEPT + document config rule |
| C3 | 1D/2D vdiff zeroing (Athena) vs always-3D (GAMER) | Not fixable sensibly | Comparing vs Athena 1D/2D runs with oblique B | ACCEPT + fix comparison protocol |
| C4 | Defaults: CR_SOURCE=0 vs src_flag=1; CR_VMAX 1e2 vs 1.0 | Trivial | Only if users rely on defaults | Document in PR |
| C5 | CR_SIGMA input units (paper sigma' vs Athena internal) | None (deliberate) | User inputs | Document prominently |
| C6 | B.C. applied once/step vs twice/step | Hard (structural) | Domain boundaries | **ACCEPT** (already agreed) |
| C7 | Double precision not enforced | Trivial | Validation runs | Add check/warning |

## Suggested order of work

### Tier 1 — do before the PR (improves the tests you already have)

1. **B1 (TINY_NUMBER)** — Define a CR-module constant equal to Athena's `1.0e-20` and use it
   for every threshold inside the CR kernel (dpc_sign deadband, Ec floor value, va/btot
   checks, HLLE degenerate check, source floor). Do **not** change GAMER's global
   `TINY_NUMBER` — that would ripple through the whole framework. This is the cheapest fix
   with the largest expected payoff: it directly targets the sgn-amplifier residual
   (triangular test ~1e-8) and makes the Ec floor value identical.
2. **C1 (fallback semantics)** — While in there, reproduce Athena's exact fallback branches:
   `DefaultStreaming` leaves sigma_adv *unchanged* when va<=TINY (GAMER sets max_opacity),
   and Athena has no `Ec > TINY` guard (safe to drop once the floor is 1e-20, so the
   1/Ec blow-up matches Athena's).
3. **B7a (ADV_* advection)** — Stop the hydro solver from advecting ADV_SIGMA/ADV_VX/VY/VZ:
   uncomment/enable the skip in the HLLE passive loop **and** make sure those flux slots are
   explicitly zeroed (or the fields skipped in the flux-divergence updates), so the stored
   fields are a well-defined pass-through instead of advected junk. Cheap, removes a
   contamination source at patch edges, and makes B7b easier to reason about.
4. **B8 guard** — Add a runtime check forbidding (or warning about) `OPT__1ST_FLUX_CORR`
   with CR_STREAMING, since the fallback re-updates CR fields with hydro-only physics.
5. **C7** — Warn (Aux_Check_Parameter) if CR_STREAMING is compiled in single precision.

After Tier 1, rerun the standard comparison suite; triangular/plateau residuals should
tighten and patch-edge noise should drop.

### Tier 2 — "coupled-parity bundle 1" (needed only for live-gas matching)

Do these together; they all live in the two source functions and only show up when the gas
moves or CR_SOURCE=1. Since I don't know which science runs are planned, the safe default
is to do them — they are small and the data is already available in the kernel:

6. **B2 (gas state)** — Read rho/momentum from the stage-updated state (`out_con` at half
   step, `g_Output` at full step) instead of `g_ConVar_In` / `g_PriVar_Half`, matching
   Athena's use of post-transport `u`. Localized swap of the read source.
7. **B3 (half-step back-reaction)** — Replicate the full-step back-reaction block in the
   half-step source (momentum + energy on `out_con`, with the 0.5*dt source), matching
   Athena's stage-agnostic `AddSourceTerms`.
8. **B4 (ec_source operator)** — Build grad(Pc) from the CR flux divergence (the source
   functions already receive `g_Flux*`; the index arithmetic already exists in
   `CR_UpdateStreaming`) instead of central differences of Ec. Moderate but mechanical.

### Tier 3 — "coupled-parity bundle 2" (evolving-B parity)

9. **B5 (half-step B level)** — Athena's stage-1 source uses t^n B; GAMER uses the CT
   half-step B. Fixing requires access to t^n cell-centered B inside `Hydro_RiemannPredict`.
   Caution: `g_PriVar_1PG` and `g_PriVar_Half_1PG` alias the same buffer, so you cannot just
   read t^n B from it mid-loop — recompute from `g_Mag_Array_In` faces or add a small buffer.
   Only worth it if evolving-B runs must match; otherwise document as a truncation-order
   difference (zero for static B).

10. **B6 (floor/heating order)** — Trivial to reorder to Athena's sequence (gas energy sees
    the *unfloored* new_ec, floor threshold 0.0, floor after energy update). Note this
    reproduces Athena's spurious gas heating when the floor triggers; GAMER's current order
    is better physics. Present both options to the boss: exact match vs documented
    improvement. Either way it is a one-screen change.

11. **B8 retry** — If exactness under solver retries is ever needed, re-derive ADV_* at the
    top of each retry iteration. Low priority: Athena has no analogous retry at all, so any
    retry behavior is already "beyond Athena".

### Not worth fixing — ask for acceptance (with rationale to present)

- **B7b (stale outermost ADV_* ghost ring)** — Structural: GAMER recomputes opacity per
  patch from ghost-filled data, and the outermost ring cannot be recomputed (no +/-1
  neighbors). A true fix means widening FLU_GHOST_SIZE or changing the ghost-fill machinery
  — framework-level surgery for a small extra-dissipation difference at patch edges. Same
  acceptance class as the B.C. cadence. Suggest: quantify it once (single-patch vs
  multi-patch run of the same problem, diff the seam cells) so the acceptance is informed.
- **C6 (B.C. once vs twice per step)** — Already accepted by the team; keep it on the list
  so the PR states it explicitly. Note the in-half-step `CR_UpdateOpacity` near domain
  boundaries is where this shows up for the CR module.
- **C2 (dt semantics)** — Re-engineering GAMER's per-criterion dt framework to Athena's
  single-CFL max(|v|+c_fast, vmax) is not worth it. Instead document the config rule:
  set `CR_CFL` = Athena's `cfl_number` and keep vmax >> gas speeds; note the codes diverge
  if gas speeds approach vmax. (Optional cheap improvement: include |v|+c_fast in the CR
  criterion, but that still uses a separate CFL knob.)
- **C3 (dimensionality)** — GAMER is inherently 3D; emulating Athena's 1D/2D vdiff zeroing
  makes no sense. Fix the *comparison protocol* instead: validate only against Athena runs
  with nx2,nx3 > 1 (or axis-aligned B in 1D, where the difference cancels).
- **C5 (CR_SIGMA units)** — Deliberate design choice (paper units in, xVm internally).
  Keep, but make it impossible to miss: PR description, wiki page, and the parameter table
  should all state "GAMER CR_SIGMA = Athena sigma / vmax".
- **C4 (defaults)** — Keep GAMER's defaults if preferred (CR_SOURCE off is safer for new
  users), just list the deviation from Athena's defaults in the PR notes.

## Verification strategy

- After Tier 1: rerun the 19-run CR_Streaming suite + classic comparisons (static gas);
  expect residual improvement, no regressions.
- Tier 2/3 cannot be validated by the current suites (they are static-gas by design).
  Before claiming coupled parity, add at least one step-aligned live-gas comparison
  (e.g. the paper shock test with CR_SOURCE=1, fixed-step dumping on both codes).
