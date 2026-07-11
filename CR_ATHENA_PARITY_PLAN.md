# CR Two-Moment vs Athena++: Parity Plan

Follow-up to the pre-PR review of the CR_STREAMING module. Each known GAMER-vs-Athena
difference is classified by **difficulty**, **impact on matching**, and a **recommendation**:
`FIX` (do it), `FIX-LATER` (only needed for live-gas/coupled parity), or
`ACCEPT` (structural GAMER-framework difference; ask for sign-off, like the B.C. cadence).

Item numbers (B1..B8, C1..C7) refer to the review report.

## Summary table

| # | Item | Difficulty | When it matters | Recommendation |
|---|------|-----------|-----------------|----------------|
| B1 | TINY_NUMBER: DBL_MIN vs 1e-20 | Easy | Edge cells only (see 2026-07-11 note) | **OPEN — deferred**; sites tagged `[B1]` in source |
| C1 | `va<=TINY` fallback + extra `Ec>TINY` guard | Easy | Edge cells (floored Ec, B~0) | **OPEN — deferred with B1**; sites tagged `[C1]` in source |
| B7a | ADV_* passive-advection fluxes (HLLE skip commented out) | Easy | Always (stored fields, patch edges) | **FIXED 2026-07-10** |
| B2 | Gas state fed to source solve (t^n / half-step vs stage-updated) | Easy–Moderate | Live gas only (CR_SOURCE=1 or moving gas) | **FIXED 2026-07-09** |
| B3 | Missing gas back-reaction at half step | Easy–Moderate | Live gas + CR_SOURCE=1 | **FIXED 2026-07-09** |
| B4 | ec_source: central-diff vs flux-divergence grad(Pc) | Moderate | Multi-D, moving gas, oblique B | **FIXED 2026-07-09** (full 9-term form; see B9) |
| B9 | grad(Pc) diagonal-only vs Athena's full 9-term flux divergence (new finding 2026-07-09) | Easy (mechanical) | Multi-D / oblique B with streaming | **FIXED 2026-07-09** (see Tier 2 item 9) |
| B5 | Half-step source uses half-step B vs Athena t^n B | Moderate (buffer aliasing) | Evolving B only | **FIXED 2026-07-10** (recompute t^n B from face-centered input) |
| B6 | Floor vs gas-energy update ordering in source | Trivial | Only when source drives Ec<0 with CR_SOURCE=1 | Decide with boss (match = trivial; GAMER's version conserves energy better) |
| B8 | MinMod retry reuses mutated sigma; 1st-order-flux-corr fallback | Easy (guard) / Moderate (retry) | Solver-failure paths only | **FIXED 2026-07-11** (guard + retry) |
| B7b | Stale ADV_* in outermost ghost ring | Hard (structural) | Patch-boundary cells, every step | **FIXED 2026-07-10, run-validated 2026-07-11** (FLU_GHOST_SIZE +1; 1D streaming now matches Athena to double precision) |
| C2 | dt criterion semantics | Easy (was rated Moderate; see section) | Only if gas speed ~ vmax | **FIXED 2026-07-10** (folded into hydro CFL solver; ACCEPT overridden) |
| C3 | 1D/2D vdiff zeroing (Athena) vs always-3D (GAMER) | Not fixable sensibly | Comparing vs Athena 1D/2D runs with oblique B | ACCEPT + fix comparison protocol |
| C4 | Defaults: CR_SOURCE=0 vs src_flag=1; CR_VMAX 1e2 vs 1.0 | Trivial | Only if users rely on defaults | Document in PR |
| C5 | CR_SIGMA input units (paper sigma' vs Athena internal) | None (deliberate) | User inputs | Document prominently |
| C6 | B.C. applied once/step vs twice/step | Hard (structural) | Domain boundaries | **ACCEPT** (already agreed) |
| C7 | Double precision not enforced | Trivial | Validation runs | **CLOSED 2026-07-11** (no action; double is a comparison-protocol choice, not a module requirement) |

## Suggested order of work

### Tier 1 — status after the 2026-07-11 decisions

**Context that changed the calculus:** the B7a+B7b fixes were run-validated on 2026-07-11 —
1D streaming now matches Athena to **double precision**. B7b (the stale ADV_* ghost ring),
not B1, was therefore the dominant cause of the old triangular ~1e-8 residual, and B1's
"largest expected payoff" (below) has already been captured.

1. **B1 (TINY_NUMBER)** — **OPEN, deferred** (decision 2026-07-11; fix deliberately not
   applied for now, difference accepted). Original idea kept for reference: define a
   CR-module constant equal to Athena's `1.0e-20` and use it for every threshold inside the
   CR kernel (dpc_sign deadband, Ec floor value, va/btot checks, HLLE degenerate check,
   source floor); do **not** change GAMER's global `TINY_NUMBER`. Remaining exposure is
   edge-case only: cells where |B·∇Pc|, va, or Ec fall between DBL_MIN and 1e-20 — none
   occur in the validated suites. All affected thresholds are tagged **`[B1]`** in
   `CPU_CR_TwoMoment.cpp` (with a summary block at the top of the file); if a future
   comparison shows noise at plateaus/extrema (roundoff-level ∇Pc), suspect these first.
2. **C1 (fallback semantics)** — **OPEN, deferred with B1** (decision 2026-07-11): keep
   GAMER's `Ec > TINY` guard and the max_opacity cap. The original fix ("drop the guard")
   was only safe with a 1e-20 floor — with the floor at DBL_MIN, Athena's 1/Ec blow-up
   would overflow, so the guard is mandatory as long as B1 stays open. Sites tagged
   **`[C1]`** in `CPU_CR_TwoMoment.cpp`. Consequence: floored-Ec / degenerate-B cells will
   not match Athena's sigma_adv (GAMER caps, Athena blows up ∝1/Ec — GAMER's is the safer
   behavior).
3. **B7a (ADV_* advection)** — Stop the hydro solver from advecting ADV_SIGMA/ADV_VX/VY/VZ:
   uncomment/enable the skip in the HLLE passive loop **and** make sure those flux slots are
   explicitly zeroed (or the fields skipped in the flux-divergence updates), so the stored
   fields are a well-defined pass-through instead of advected junk. Cheap, removes a
   contamination source at patch edges, and makes B7b easier to reason about.
   **DONE 2026-07-10.** Implemented by ZEROING the ADV_* flux slots inside
   `CR_TwoMomentFlux_HalfStep/FullStep()` (new step 9 in both), NOT by enabling the HLLE
   skip: the zeroing is Riemann-solver-agnostic (the suites use HLLD/HLLE), and skipping
   inside a solver would leave uninitialized flux slots that the `RSOLVER_RESCUE` NaN scan
   reads. The CR flux loops cover exactly the faces read by `Hydro_RiemannPredict()` /
   `Hydro_FullStepUpdate()`, so dflux[ADV_*] = 0 everywhere → exact pass-through for any
   solver. The dead commented-out skip in `CPU_Shared_RiemannSolver_HLLE.cpp` is removed.
   Stored ADV_* now hold well-defined values (the half-step `CR_UpdateStreaming` output at
   PS2 cells) but remain diagnostic-only after B7b — they will NOT match Athena's dumped
   end-of-step `DefaultOpacity` sigma, so compare Ec/Fc/gas across codes, not sigma.
4. **B8 (guard + retry)** — **FIXED 2026-07-11**; see the dedicated section below.
5. **C7** — **CLOSED 2026-07-11, no action.** The CR module does not require double
   precision by nature; `--double=true` was only ever a comparison-protocol choice for the
   Athena validation runs. Too basic to warrant even a user-facing note.

### B8 fix (2026-07-11) — solver-failure paths (compile-verified, per decision: no runs)

Finding 8 had two INDEPENDENT halves; fixing one does not cover the other:

- **Guard (host-side 1st-order flux correction):** `Flu_Close()`'s fallback recomputes
  failing cells with pure-hydro 1st-order Riemann fluxes — CR_E/CR_F* become passively
  advected scalars (no two-moment flux, no source) and ADV_* gets re-advected (also voiding
  the B7a pass-through there). Critically, `OPT__1ST_FLUX_CORR` **defaults ON** for MHD
  (`FIRST_FLUX_CORR_3D`, `Init_ResetParameter.cpp`), so every CR_STREAMING run had it armed.
  Fix follows the existing SRHD precedent: `Init_ResetParameter.cpp` now defaults
  `OPT__1ST_FLUX_CORR = FIRST_FLUX_CORR_NONE` under CR_STREAMING when the user left it at
  -1; if the user enables it explicitly, `Aux_Check_Parameter.cpp` fires a warning stating
  the cost (CR fields updated with hydro-only physics wherever the correction triggers).
- **Retry (in-kernel MinMod loop, `MINMOD_MAX_ITER > 0` only; default 0 = off):** the retry
  re-runs reconstruction→fluxes→streaming→source with full CR physics, but the previous
  iteration's `CR_UpdateStreaming()` had overwritten ADV_* in `g_PriVar_Half` in-place, so
  retries saw DefaultStreaming values where the first attempt saw DefaultOpacity values.
  Fix: at the top of the retry body (`CPU_FluidSolver_MHM.cpp`, MHM_RP do-loop), guarded by
  `Iteration > 0`, re-issue the exact end-of-`Hydro_RiemannPredict()` `CR_UpdateOpacity()`
  call. `CR_UpdateOpacity()` reads only DENS/CR_E/B — none written inside the loop — so the
  restore is bitwise-exact, and the `Iteration > 0` guard makes the non-retry path provably
  unchanged (no revalidation of existing suites needed). GPU via the `.cu` symlink
  (`__syncthreads()` after the call; `Iteration` is block-uniform so the branch is safe);
  no new memory, no signature changes. Note Athena has no retry at all, so this fix buys
  GAMER *iteration-independence* (self-consistency), not Athena parity — for strict parity
  keep `MINMOD_MAX_ITER=0` (the default).

Compile-verified (cpu_CR_Streaming, other_tests/gpu_CR_Streaming_x, cpu_CR_Classic_Diffusion
— the last confirms the non-CR_STREAMING paths of `Init_ResetParameter`/`Aux_Check_Parameter`
still build). Per decision, no runtime tests: the guard is init-time logic and the retry is
bitwise-inert unless a full-step failure occurs with `MINMOD_MAX_ITER > 0`.

### B7b fix (2026-07-10) — stale ADV_* ghost ring eliminated by widening FLU_GHOST_SIZE

Originally classified **ACCEPT** (structural); overridden by decision 2026-07-09 — exact
patch-size independence is wanted (result with PATCH_SIZE=8 must equal PATCH_SIZE=16).

Mechanism recap: the start-of-step `CR_UpdateOpacity()` recomputes ADV_* only on FLU_NXT
layers [1, FLU_NXT-2]; layer 0 keeps ghost-filled *stored* values (a neighbor's end-of-step
ADV_*, not what a wider patch would recompute there). The contamination penetrates exactly
3 layers per step — ring(0) → half-step cell(1) via vdiff at face 0 → {ADV layer 2 via the
in-`Hydro_RiemannPredict` opacity central diff; PLM slope of the outermost FC cell} →
full-step flux at the PS2-edge face → PS2 edge cell — and the old MHM_RP ghost size is
exactly 3 (`2 + LR_GHOST_SIZE`), so it reached the output every step.

**Fix: `FLU_GHOST_SIZE = 3 + LR_GHOST_SIZE` under CR_STREAMING** (one extra layer; the same
+1 also suffices for PPM, whose deeper stencil is offset by the deeper FC region). Every
contamination path now dies one layer short of the PS2 output.

Implementation (all sites; compile-verified 2026-07-10 with the PLM config, the PPM+HLLE
config, and a classic-CR config — the last confirms non-CR_STREAMING builds keep ghost
size 3 and are bit-identical):
- `include/Macro.h`: MHM_RP branch gains an `#ifdef CR_STREAMING` giving `3 + LR_GHOST_SIZE`
  (verified: FLU_GHOST_SIZE=4, FLU_NXT=24 at PS1=8 with PLM).
- `CPU_Shared_DataReconstruction.cpp` (both variants) and `CPU_FluidSolver_MHM.cpp`
  (`OffsetPri` for `MHD_ComputeElectric`): the FC↔half-step centering offsets hardcoded
  `LR_GHOST_SIZE`; now computed generically — `(NIn - N_FC_VAR)/2` and
  `(N_HF_VAR - N_FC_VAR)/2` — value-identical for every non-CR build (MHM/MHM_RP/CTU).
- `Aux_Check_Parameter.cpp`: compile-time assert `FLU_GHOST_SIZE == 3 + LR_GHOST_SIZE`
  under CR_STREAMING; `OPT__LR_LIMITER = LR_LIMITER_EXTPRE` is now an ERROR with
  CR_STREAMING (its ±2 stencil reaches the stale ring again and would void the guarantee);
  the "reduce FLU_GHOST_SIZE for performance" warnings are suppressed for CR_STREAMING.
- No CR-kernel changes needed: all CR offsets were already macro-generic
  (`fc2pvar_offset`, `pvar_offset`, `half_offset`), and `N_FC_VAR`/`N_FL_FLUX`/PS2 mappings
  are unchanged — only the half-step region (`N_HF_VAR = FLU_NXT-2`) widens automatically.

Consequences:
- Per-patch ADV_* recomputation is now equivalent to a global one: **uniform-grid results
  must be bitwise identical across patch sizes** (PS 8 vs 16 — the acceptance test; dt is a
  min-reduction, hence order-exact). AMR coarse–fine boundaries still interpolate ghost
  data — inherent to AMR, not B7b.
- Stored ADV_* fields no longer influence the evolution AT ALL (they only ever entered
  through the ring); they are diagnostic-only (dumps/restarts).
- Cost: FLU_NXT 22→24 at PS1=8 (~+30% ghost-fill/half-step work, +33% MPI halo width);
  the full step is unchanged; applies to CR_STREAMING builds only.
- Archived GAMER-vs-Athena results will shift at patch seams (in the improving direction);
  rerun before comparing.

### C2 fix (2026-07-10) — Athena dt semantics folded into the hydro CFL solver

Originally classified **ACCEPT** ("re-engineering GAMER's per-criterion dt framework is not
worth it"); overridden — the Moderate rating collapses to Easy on one observation: **CR_VMAX
is spatially constant, so `max(vmax, ·)` commutes with the patch-wide reduction**. Athena's
per-cell criterion `min_cells dh/max(vmax, |v|+c_fast)` equals `dh/max(vmax, MaxCFL)`, and
`MaxCFL` (the patch max of |v|+c_fast) is already computed and reduced by the fluid dt
solver. No per-cell loop, no new reduction, no framework change.

Implementation (run-validated, not just compile-verified):
- `CPU_dtSolver_HydroCFL.cpp` (GPU covered via the `.cu` symlink): under CR_STREAMING the
  per-patch output becomes
  `g_dt_Array[p] = FMIN( dhSafety/MaxCFL, MicroPhy.CR_cfl*dh / FMAX(MicroPhy.CR_vmax, MaxCFL) )`
  — the same pattern CR_DIFFUSION already uses to fold a second criterion into this solver.
  GPU impact: ~2 flops in the existing thread-0 epilogue; no new syncs, memory, or signature
  changes.
- `Mis_GetTimeStep.cpp` criterion 1.10 (flat `CR_CFL*dh/CR_VMAX`) is **kept, unchanged in
  value**, as a deliberate fallback: `OPT__FREEZE_FLUID` resets criterion ONE to HUGE_NUMBER,
  which would otherwise leave frozen-fluid CR runs with no dt limit. It is redundant
  (always >= the folded term) whenever the fluid criterion is active.
- `Aux_Check_Parameter.cpp` warns when `CR_CFL > DT__FLUID`.

**New config rule (supersedes "keep vmax >> gas speeds"):** with `CR_CFL <= DT__FLUID`, the
folded term always wins the min, so GAMER's dt == Athena's `cfl_number*dh/max(|v|+c_fast, vmax)`
with `CR_CFL` playing the role of Athena's `cfl_number` — at **all** gas speeds.

Remaining (accepted) dt caveats:
- Athena's dt fast speed augments the normal component: `bx = bcc + |b_face - bcc|`
  (`hydro/new_blockdt.cpp`); GAMER uses the pure cell-centered B. Identical for B uniform
  within a cell; otherwise a pre-existing *global hydro-CFL* difference (any MHD comparison),
  and it only matters where the gas side wins the max.
- dt facet of C3: Athena 1D/2D runs skip the collapsed directions in the dt min; GAMER is
  always 3D. Same acceptance/protocol as C3.
- At Step 0 the fluid part uses `DT__FLUID_INIT`; keep it >= CR_CFL (or accept a smaller
  first step vs Athena).
- `Record__TimeStep` column semantics: the "Hydro_CFL" column now includes the CR cap under
  CR_STREAMING (it ties with "CR_Stream" whenever vmax dominates); "CR_Stream" remains the
  flat light-crossing value.

Verification (2026-07-10, runs in `bin/agents/cr_c2_dt_semantics/`, triangular test, 10 steps):
- **No-op check** (CR_VMAX=1e2 >> gas fast speed 1.633): selected dt bitwise unchanged
  (2.3437500e-5 every step, time history identical); only the Hydro_CFL display column
  changed (3.8273277e-3 -> 2.3437500e-5, the tie). Confirms all archived/standard suites
  are unaffected.
- **Active-regime check** (CR_VMAX=0.5 < 1.633): post-fix dt = 1.4352479e-3 =
  `0.3*dh/1.6329932` — Athena's formula with cfl_number=0.3, matching the analytic
  perpendicular fast speed sqrt(8/3) to all printed digits; pre-fix gave 3.8273277e-3
  (old mixed semantics).
- Compile checks: cpu_CR_Streaming, other_tests/gpu_CR_Streaming_x (GPU), and
  cpu_CR_Classic_Diffusion (non-CR_STREAMING path untouched).

### Tier 2 — "coupled-parity bundle 1" (needed only for live-gas matching)

**STATUS: DONE 2026-07-09** (all three, in `CPU_CR_TwoMoment.cpp` only; no signature changes,
GPU covered via the `.cu` symlink; compile-verified with the cpu_CR_Streaming config; not yet
validated by a live-gas comparison run — see Verification strategy).

Do these together; they all live in the two source functions and only show up when the gas
moves or CR_SOURCE=1. Since I don't know which science runs are planned, the safe default
is to do them — they are small and the data is already available in the kernel:

6. **B2 (gas state)** — Read rho/momentum from the stage-updated state (`out_con` at half
   step, `g_Output` at full step) instead of `g_ConVar_In` / `g_PriVar_Half`, matching
   Athena's use of post-transport `u`. Localized swap of the read source.
   **DONE 2026-07-09.** Implementation notes: rho is guarded with `FMAX(rho, TINY_NUMBER)`
   because both stage-updated states are density-unfloored at the point of the source call
   (half step: floor applied after the source in `Hydro_RiemannPredict()`; full step:
   `Hydro_FullStepUpdate()` defers flooring to `Flu_Close()`). Athena instead applies its EoS
   density floor *inside* the source (`cr_source.cpp: rho = max(rho, rho_floor)`), so cells
   where either floor triggers will not match — this ties into the still-open B6 decision.
   The ec_source velocities deliberately KEEP the pre-stage state (`g_ConVar_In` at half step,
   `g_PriVar_Half` at full step): Athena evaluates `ec_source_` in `CalculateFluxes()` from
   the pre-stage `w` while the implicit solve uses stage-updated `u`.
7. **B3 (half-step back-reaction)** — Replicate the full-step back-reaction block in the
   half-step source (momentum + energy on `out_con`, with the 0.5*dt source), matching
   Athena's stage-agnostic `AddSourceTerms`.
   **DONE 2026-07-09.** Mirrors GAMER's full-step ordering (floor new_ec first, then kick the
   gas), which deviates from Athena's order — that ordering question remains item B6. The
   kicked half-step gas then feeds the MinEint floor, Con2Pri, and the full-step
   reconstruction, matching Athena's stage-1 flow.
8. **B4 (ec_source operator)** — Build grad(Pc) from the CR flux divergence (the source
   functions already receive `g_Flux*`; the index arithmetic already exists in
   `CR_UpdateStreaming`) instead of central differences of Ec. Moderate but mechanical.
   **DONE 2026-07-09**, using the full 9-term discretization (see item 9 / B9). Half step
   uses the previously-unused `idx_flux`/`didx_flux` arguments on `g_Flux_Half`; full step
   computes the flux index internally (`i_out+1` mapping, identical to
   `Hydro_FullStepUpdate()` with MHD).
9. **B9 (grad(Pc): full 9-term flux divergence)** — New finding during the Tier-2 work
   (2026-07-09): Athena's `grad_pc_(n)` is the FULL divergence of the flux vector of the
   CRF(n+1) equation — it sums the flux differences of ALL THREE directions per component
   (9 terms total; note the `+=` in the nx2>1/nx3>1 blocks of `cr_transport.cpp:366-403`),
   and the transverse fluxes are generically nonzero through the HLLE upwind/dissipation
   terms (`-bm*Fc_L`, `-bp*Fc_R`) even though the off-diagonal Eddington factors vanish.
   Athena feeds the SAME `grad_pc_` array to both `DefaultStreaming` (sigma_adv/v_adv via
   `b_grad_pc = B·grad_pc`) and `ec_source_`, so both consumers need all 9 terms. GAMER's
   `CR_UpdateStreaming` had kept only the 3 diagonal terms — zero effect for the 1D
   axis-aligned suites (transverse fluxes uniform), but breaks machine-precision parity in
   multi-D / oblique-B streaming (possibly part of the circle test's 2.9e-4 residual).
   **DONE 2026-07-09** (initially deferred pending verification of the Athena source, then
   confirmed — the x2/x3 accumulations are easy to misread — and fixed the same day):
   `CR_UpdateStreaming` and BOTH ec_source blocks now compute the full 9-term divergence.
   All required transverse fluxes were verified to be computed and in range at every call
   site (no ghost widening, no new GPU syncs).

### Tier 3 — "coupled-parity bundle 2" (evolving-B parity)

9. **B5 (half-step B level)** — Athena's stage-1 source uses t^n B; GAMER uses the CT
   half-step B. Fixing requires access to t^n cell-centered B inside `Hydro_RiemannPredict`.
   Caution: `g_PriVar_1PG` and `g_PriVar_Half_1PG` alias the same buffer, so you cannot just
   read t^n B from it mid-loop — recompute from `g_Mag_Array_In` faces or add a small buffer.
   Only worth it if evolving-B runs must match; otherwise document as a truncation-order
   difference (zero for static B).
   **FIXED 2026-07-10** (compile-verified, not run-validated) via the "recompute from
   `g_Mag_Array_In`" option: `Hydro_RiemannPredict()` gains a `g_FC_B_In` parameter (the t^n
   face-centered B, `g_Mag_Array_In[P]`) and recomputes the t^n cell-centered B per cell with
   `MHD_GetCellCenteredBField()` — bitwise identical to the t^n B used by the half-step
   fluxes (same array, same averaging, same indices). `CR_TwoMomentSource_HalfStep()` takes
   the result as a new `B_n[]` parameter and uses it for the rotation angles instead of
   `OneCell[MAG_OFFSET+*]`; a single set of angles feeds the implicit solve AND both
   ec_source rotations, so the one switched read fixes all B consumers at once (matching
   Athena, where stage 1 uses t^n `b_angle` everywhere). The full-step source keeps the
   half-step B (= Athena's stage-2 bcc, already correct). The buffer-aliasing risk was fully
   sidestepped: `g_Mag_Array_In` is a read-only (`const`) kernel argument that nothing in the
   kernel writes, so the fix needs no new global/shared memory, no extra `__syncthreads()`,
   and no kernel-signature/CUAPI changes (~6 extra global loads per cell). Compile-verified
   with cpu_CR_Streaming, cpu_CR_Classic_Diffusion (non-CR_STREAMING path, signature change
   only), and other_tests/gpu_CR_Streaming_x (GPU + COSMIC_RAY + CR_STREAMING combined).

10. **B6 (floor/heating order)** — Trivial to reorder to Athena's sequence (gas energy sees
    the *unfloored* new_ec, floor threshold 0.0, floor after energy update). Note this
    reproduces Athena's spurious gas heating when the floor triggers; GAMER's current order
    is better physics. Present both options to the boss: exact match vs documented
    improvement. Either way it is a one-screen change.

11. **B8 retry** — **FIXED 2026-07-11** exactly as described here (re-derive ADV_* at the
    top of each retry iteration); see the dedicated B8 section under Tier 1. Athena has no
    analogous retry at all, so this is a self-consistency fix, "beyond Athena".

### Not worth fixing — ask for acceptance (with rationale to present)

- ~~**B7b (stale outermost ADV_* ghost ring)**~~ — was listed here (same acceptance class
  as the B.C. cadence), but the ACCEPT was overridden and it is **FIXED 2026-07-10** by
  widening FLU_GHOST_SIZE; see the dedicated section above.
- **C6 (B.C. once vs twice per step)** — Already accepted by the team; keep it on the list
  so the PR states it explicitly. Note the in-half-step `CR_UpdateOpacity` near domain
  boundaries is where this shows up for the CR module.
- ~~**C2 (dt semantics)**~~ — was listed here ("re-engineering the per-criterion dt
  framework is not worth it"), but no re-engineering was needed: since CR_VMAX is spatially
  constant, Athena's per-cell max folds into the already-reduced MaxCFL of the fluid dt
  solver. **FIXED 2026-07-10**; see the dedicated section above. Config rule is now
  `CR_CFL <= DT__FLUID` (with CR_CFL = Athena's `cfl_number`) for exact dt parity at all
  gas speeds.
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
- Tier-2 landing note (2026-07-09): Tier 2 (+ B9) was implemented BEFORE Tier 1 (user
  request) and is compile-verified only. For exactly static gas (v=0, zero mass/momentum
  fluxes), B2/B3/B4 are mathematical no-ops: B2's rho/v feed only v-dependent terms, B3 is
  gated by CR_SOURCE, and B4's new gradient is multiplied by v_perp. B9 is NOT a static
  no-op in multi-D: the 9-term grad_pc changes sigma_adv/v_adv wherever transverse CR fluxes
  vary (e.g. the circle test — residuals expected to shift, presumably toward Athena);
  1D axis-aligned suites are unaffected. Runs with live gas or CR_SOURCE=1 (e.g. the
  paper-shock 4.2.2 reproduction) WILL shift at O(dt) — rerun before comparing against
  archived results.
- B7a/B7b landing note (2026-07-10): compile-verified only (PLM, PPM+HLLE, and classic-CR
  configs). **RUN-VALIDATED 2026-07-11: 1D streaming now matches Athena to double
  precision** — confirming B7b was the dominant cause of the old triangular ~1e-8 residual
  (the earlier "sgn-amplifier, unfixable" attribution is superseded).
  B7b acceptance test: identical uniform-grid runs with PATCH_SIZE=8 vs 16
  (fixed-step dumps, `--double=true`) must be **bitwise identical** after the fix — and
  measurably different before it (that pre-fix diff quantifies the old seam error).
  Both fixes change results wherever patch boundaries exist, even in 1D static tests
  (B7a removes advected-junk feedback through the ring, B7b removes the ring's influence
  entirely), so patch-seam residuals vs Athena should DROP; rerun any archived comparison
  before reuse. Dumped ADV_* fields are now pass-through diagnostics and will not match
  Athena's dumped sigma — compare Ec/Fc/gas fields only.
