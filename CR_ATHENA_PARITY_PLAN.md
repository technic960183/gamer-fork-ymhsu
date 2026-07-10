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
| B7a | ADV_* passive-advection fluxes (HLLE skip commented out) | Easy | Always (stored fields, patch edges) | **FIXED 2026-07-10** |
| B2 | Gas state fed to source solve (t^n / half-step vs stage-updated) | Easy–Moderate | Live gas only (CR_SOURCE=1 or moving gas) | **FIXED 2026-07-09** |
| B3 | Missing gas back-reaction at half step | Easy–Moderate | Live gas + CR_SOURCE=1 | **FIXED 2026-07-09** |
| B4 | ec_source: central-diff vs flux-divergence grad(Pc) | Moderate | Multi-D, moving gas, oblique B | **FIXED 2026-07-09** (full 9-term form; see B9) |
| B9 | grad(Pc) diagonal-only vs Athena's full 9-term flux divergence (new finding 2026-07-09) | Easy (mechanical) | Multi-D / oblique B with streaming | **FIXED 2026-07-09** (see Tier 2 item 9) |
| B5 | Half-step source uses half-step B vs Athena t^n B | Moderate (buffer aliasing) | Evolving B only | **FIXED 2026-07-10** (recompute t^n B from face-centered input) |
| B6 | Floor vs gas-energy update ordering in source | Trivial | Only when source drives Ec<0 with CR_SOURCE=1 | Decide with boss (match = trivial; GAMER's version conserves energy better) |
| B8 | MinMod retry reuses mutated sigma; 1st-order-flux-corr fallback | Easy (guard) / Moderate (retry) | Solver-failure paths only | Guard now, document retry |
| B7b | Stale ADV_* in outermost ghost ring | Hard (structural) | Patch-boundary cells, every step | **FIXED 2026-07-10** (FLU_GHOST_SIZE +1; ACCEPT overridden) |
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
4. **B8 guard** — Add a runtime check forbidding (or warning about) `OPT__1ST_FLUX_CORR`
   with CR_STREAMING, since the fallback re-updates CR fields with hydro-only physics.
5. **C7** — Warn (Aux_Check_Parameter) if CR_STREAMING is compiled in single precision.

After Tier 1, rerun the standard comparison suite; triangular/plateau residuals should
tighten and patch-edge noise should drop.

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

11. **B8 retry** — If exactness under solver retries is ever needed, re-derive ADV_* at the
    top of each retry iteration. Low priority: Athena has no analogous retry at all, so any
    retry behavior is already "beyond Athena".

### Not worth fixing — ask for acceptance (with rationale to present)

- ~~**B7b (stale outermost ADV_* ghost ring)**~~ — was listed here (same acceptance class
  as the B.C. cadence), but the ACCEPT was overridden and it is **FIXED 2026-07-10** by
  widening FLU_GHOST_SIZE; see the dedicated section above.
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
  configs). B7b acceptance test: identical uniform-grid runs with PATCH_SIZE=8 vs 16
  (fixed-step dumps, `--double=true`) must be **bitwise identical** after the fix — and
  measurably different before it (that pre-fix diff quantifies the old seam error).
  Both fixes change results wherever patch boundaries exist, even in 1D static tests
  (B7a removes advected-junk feedback through the ring, B7b removes the ring's influence
  entirely), so patch-seam residuals vs Athena should DROP; rerun any archived comparison
  before reuse. Dumped ADV_* fields are now pass-through diagnostics and will not match
  Athena's dumped sigma — compare Ec/Fc/gas fields only.
