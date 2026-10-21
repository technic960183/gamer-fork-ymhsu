# cr_streaming branch: verbose / stale comment review (2026-10-05)

Scope: every comment added on `cr_streaming` since merge-base ee11f8a8a. Line numbers are HEAD 58db8b061.
Tags: RESTATE (narrates code) · FILLER (empty words) · HISTORY (dev story / Athena-parity log) ·
OVERJUSTIFY (paragraph where a clause will do) · PARAM (bloated param entry) · DECOR (CAPS, NOTE:, unicode) ·
SPECULATIVE · STALE (comment is WRONG vs code — fix these first).
"(yours)" = came from the user's own commit 3d46112f8.

---------------------------------------------------------------------------------------------------
## 0. Fix first — comments that are WRONG (STALE)

| where | problem |
|---|---|
| configure.py L707-710 | help says it can be built with `--cosmic_ray`/`--eos=COSMIC_RAY`; Aux_Check_Parameter L1956 `#error`s unless EOS_GAMMA (constraint on L706 also still allows COSMIC_RAY) |
| Mis_GetTimeStep.cpp L241-244, dtSolver L207-215 | "kept as fallback for OPT__FREEZE_FLUID" is false: InvokeSolver.cpp L104 sets fluid dt=0 → CR frozen too. Criterion Ten is fully redundant |
| Aux_Check_Parameter.cpp L1984-1985 (yours) | missing space → prints "will not begood" |
| CR_TwoMoment L297-298 + MHM L1095 | "used in the NEXT timestep's flux" — wrong; MHM L492 feeds this step's half-step flux, L1107 feeds this step's full-step flux |
| CR_TwoMoment L131 | "sigma_adv[0] parallel, [1,2] perpendicular" — sigma_adv is a scalar now |
| CR_TwoMoment L132 | "Invoked by CR_UpdateStreaming()" — also CR_UpdateOpacity() since 3715eb58d |
| CR_TwoMoment L201, L307 | stride list includes PS2 — no caller passes it |
| CR_TwoMoment L606-607, L745 | `g_Con_Var` (real name `g_ConVar`); "hydrodynamic fluxes for adding the cosmic-ray diffusive fluxes" — copied from CR_AddDiffuseFlux; this overwrites two-moment fluxes |
| CR_TwoMoment L684-685, L842-843 | "computed at init or previous flux call" — set by CR_UpdateOpacity() just before |
| CR_TwoMoment L704, L867 | "add fluxes" — code overwrites (L867 even says "replace, not add") |
| CR_TwoMoment L938-939 | "loop above ... out_con" — leftover from when this was inline in the caller |
| CR_TwoMoment L980-981, L1244-1245 | "updated by flux function" — written by CR_UpdateStreaming(); "Python/Athena++ pattern" refers to a prototype |
| CR_TwoMoment L1015, L1276 | "perpendicular components have no streaming contribution" — code zeroes gas velocity too; real reason: perp work goes via ec_source |
| CR_TwoMoment L500, L504, L530 | "CR flux" (it's reduced Fc/vmax); "0=CRF1.." (means x/y/z); "clamp to light speed" (it's vmax/sqrt(3)) |
| CR_TwoMoment L412 | Athena "lines 77-170" — points into a locally modified Athena tree |
| MHM L702-704 | "step 1-a-1" — that call is unlabeled, between 1-a-1 and 1-a-2 |
| Microphysics.h L26-27 (yours) | stray ")", undefined "Vm", two lines worded differently |
| Init_Field.cpp L152 | "CR flux fields" — block also adds CR_E |
| Main.cpp L373 | "cosmic-ray streaming" → "cosmic-ray two-moment" |
| TestProb L5-10 | "each example directory carries..." — only one example dir exists |
| TestProb L20, Input__TestProb L11 | "frozen gas" — gas is live when CR_SOURCE=1 |
| TestProb L32, Input__TestProb L18 | Dir doc omits 3 = (1,1,1) |
| TestProb L401-410 | "see the example input" — no example input for mode 7 |
| TestProb L465-475 | ignores FcInit=0; shipped Input__Parameter has CR_SIGMA=1e6, not 10 |
| TestProb L500-501 | omits /CR_VMAX and FcInit=0 |
| TestProb L733-735 | "gas keeps standard outflow" — code writes the IC |
| TestProb L118, L829 | template instruction "replace HYDRO by the target model" — upstream tests delete it |
| CR_Diffusion L237-239 | references non-existent example/test_problem/Hydro/CR_Classic_Diffusion |
| Input__Parameter L138/139/149 | `[0/1]` where the default goes (defaults: 0 / 1 / 0) |
| Input__Parameter L142 | missing default [1e2] |

---------------------------------------------------------------------------------------------------
## 1. src/Model_Hydro/CPU_Hydro/CPU_CR_TwoMoment.cpp

### InvRotateVec / CR_ComputeBFieldAngles
- L56 RESTATE `1. Applies R1^-1 * R2^-1` → DELETE (L78 + matrices say it)
- L108-109, L116-117 RESTATE `sin(theta) = |Bxy|/|B|` ... → DELETE (code + RotateVec params say it)

### CR_UpdateStreaming_OneCell
- L131 STALE → `1. sigma_adv is the opacity parallel to B`
- L132 STALE → `2. Invoked by CR_UpdateStreaming() and CR_UpdateOpacity()`

### CR_UpdateStreaming
- L189 RESTATE `Loops over interior cells and calls CR_UpdateStreaming_OneCell() for each` → `1. Invoke CR_UpdateStreaming_OneCell()`
- L190 FILLER `Works for both half-step and full-step by using appropriate parameters` → `2. Work for both the half and full steps` (or DELETE — NFlux/NVar_* entries already say it)
- L191 FILLER `Invoked only after flux calculation with grad_pc computed from flux divergence` → `3. grad(Pc) = div(CR_F* fluxes)/vmax --> invoke after the CR flux calculation`
- L192-194 SPECULATIVE `Both callers currently pass the same array ... kept for a future dedicated ADV_* buffer` → DELETE (all call sites pass NVar_Out==NVar_In==NVar_B, out_offset==in_offset; either simplify the API or drop the note)
- L196 PARAM `Array to store the updated opacity (ADV_SIGMA, ADV_VX, ...)` → `Array to store the output ADV_*`
- L197 PARAM → `Array storing the input DENS and CR_E`
- L198, L306 PARAM `... [3][ CUBE(N) ]` → `Array storing the input cell-centered B field` (N undefined)
- L199 PARAM → `Array storing the input face-centered fluxes`
- L201 STALE → `(FLU_NXT or N_HF_VAR)`
- L209 PARAM `Modified ADV_SIGMA, ADV_VX, ADV_VY, ADV_VZ in g_Output` → `g_Output[ADV_*]`
- L242, L248 RESTATE `(with out_offset)` / `(with in_offset)` → drop parentheticals
- L254 RESTATE `(may differ from input index if NVar_B != NVar_In)` → drop

### CR_UpdateOpacity
- L290-291 DECOR `... using CENTRAL DIFFERENCE gradient of CR pressure` → `Update sigma_adv and v_adv using the central-difference grad(Pc)`
- L293 RESTATE `Based on Athena++ DefaultOpacity()` → DELETE (Reference L318 says it)
- L294 DECOR `This function uses central differences ..., NOT flux divergence` → DELETE (repeats Description)
- L295-296 OVERJUSTIFY `Should be called AFTER the source term (like Athena's ...)` → `Invoke after the CR source term, as Athena++ DefaultOpacity()`
- L297-298 STALE `used in the NEXT timestep's flux calculation` → DELETE or `Output is used by the next CR flux calculation`
- L299-302 OVERJUSTIFY conserved/primitive paragraph → `g_Output[] may be conserved or primitive: only DENS, CR_E, and ADV_* are accessed` (MHM L1099-1101 repeats it — delete there)
- L304-305 PARAM → `Array storing the input DENS/CR_E and the output ADV_*`
- L307 STALE `Size for computing output indices (..., or PS2)` → `Stride of the output index (FLU_NXT or N_HF_VAR)`
- L308 PARAM `Size for computing input indices for neighbor access` → `Stride of the input index (FLU_NXT or N_HF_VAR)`
- L316 PARAM → `g_Output[ADV_*]`
- L362-363 RESTATE (optional) → `grad(Pc) by central difference with Pc = Ec/3`

### CR_ComputeVdiff
- L394 DECOR `at a CELL CENTER` → `Compute the cell-centered diffusion velocity v_diff along a given direction`
- L396 RESTATE `Based on Athena++ cr_transport.cpp CalculateFluxes()` → DELETE (move to Reference)
- L397 RESTATE `Computes v_diff in B-aligned frame, then rotates` → DELETE (L470 says it)
- L398 RESTATE `Returns v_diff component along ...` → DELETE (repeats Return)
- L399 RESTATE `This should be called separately for left and right cells` → DELETE
- L412 STALE `cr_transport.cpp lines 77-170` → `cr_transport.cpp CRIntegrator::CalculateFluxes()`
- L430-431 RESTATE → `opacities in the B-aligned frame (x' // B)`

### CR_ComputeHLLEFlux
- L496 → move to `Reference : Athena++ src/cr/integrators/cr_flux.cpp CRIntegrator::CRFlux()`
- L497 RESTATE `fdir is the flux direction (0=x, 1=y, 2=z)` → DELETE (repeats param)
- L500 STALE `CR flux (left/right)` → `reduced CR flux Fc/vmax (left/right)`
- L504 STALE `(0=CRF1, 1=CRF2, 2=CRF3)` → `(0/1/2 = x/y/z)`
- L530 STALE `clamp to light speed` → `clamp to vmax*sqrt(edd)`

### CR_TwoMomentFlux_HalfStep
- L597-598 OVERJUSTIFY `Donor-cell reconstruction matches Athena++'s forced first-order ...` → `1. First order, as Athena++ VL2 stage 1`
- L600-602 OVERJUSTIFY `Overwrites the CR slots ... (sigma_adv/v_adv are work arrays ... --> matches Athena++, ...)` → `3. Overwrite the CR_* fluxes and zero the ADV_* fluxes (ADV_* are work arrays, not advected)`
- L606 STALE `g_Con_Var` → `g_ConVar`
- L607 STALE → `Array storing the face-centered fluxes (CR_* and ADV_* overwritten)`
- L663 RESTATE `get CR variables (left and right states) - using donor cell (1st order)` → `1. get CR variables`
- L684-685 STALE+DECOR `READ sigma_adv from stored fields (computed at init or previous flux call) / Use cell-centered values directly (not averaged)` → `4. get sigma_adv (set by CR_UpdateOpacity())`
- L689 RESTATE `5. get CR parameters` → DELETE, renumber
- L692-695 HISTORY `compute v_diff at cell centers using the new helper function / vdiff_L is computed from LEFT cell's ... / This matches Athena++ pattern ...` → `6. get v_diff of the left/right cells`
- L704 STALE `add fluxes to output arrays` → `8. overwrite the CR fluxes`
- L710-714 OVERJUSTIFY 5-line ADV_* zeroing rationale → `9. zero the ADV_* fluxes so ADV_* are not advected`
- (non-comment) trailing whitespace L415, L679

### CR_TwoMomentFlux_FullStep
- L733-734 RESTATE/STALE `Uses g_FC_Var for reconstructed ... / Uses g_PriVar_Half for ... (B, sigma_adv)` → `1. vdiff uses the cell-centered g_PriVar_Half[], not the reconstructed states`
- L735 RESTATE `Uses NFlux as the stride` → DELETE
- L736-738 OVERJUSTIFY (same as half-step L600-602) → `2. Overwrite the CR fluxes and zero the ADV_* fluxes (work arrays, not advected)`
- L744 PARAM `(for B field and sigma_adv used by vdiff computation)` → DELETE; L743 → `Array storing the input cell-centered, half-step primitive variables`
- L745 STALE → `Array storing the full-step fluxes (CR and ADV_* slots are overwritten)`
- L765-766 RESTATE `g_Var has size N_FC_VAR ... N_FC_VAR = PS2+2 = 18` → DELETE (also hardcodes PS=8)
- L800, L806, L812 RESTATE → `face-centered variable index` / merge into `left/right cell-centered indices in g_PriVar_Half[]`
- L816-818 DECOR/RESTATE `from RECONSTRUCTED ... Left state: right face of left cell -> g_FC_Var[faceR][v][idx_fc] ...` → `1. get the reconstructed CR variables` + `--> L/R state = faceR of the left cell / faceL of the right cell`
- L824-825 RESTATE/DECOR → `2. get the reconstructed density and velocity`
- L831-833 OVERJUSTIFY → `3. get the cell-centered B field for vdiff`
- L842-843 STALE/DECOR → `4. get sigma_adv (set by CR_UpdateOpacity() in Hydro_RiemannPredict())`
- L850-854 HISTORY/RESTATE `using the new helper function / LEFT/RIGHT ... / This matches Athena++ pattern ...` → `6. compute vdiff from the left/right cell-centered states` (+ optional one-line Athena index map)
- L867 STALE `add fluxes ... (replace, not add, ...)` → `8. overwrite the CR fluxes`
- L873-877 OVERJUSTIFY → `9. zero the ADV_* fluxes so that Hydro_FullStepUpdate() does not advect them`

### CR_TwoMomentSource_HalfStep
- L896 RESTATE → `1. Invoked by Hydro_RiemannPredict() after the flux-divergence update`
- L897 FILLER `applies the implicit source term solve at half-step with dt_source = 0.5*dt` → DELETE
- L898 OVERJUSTIFY/DECOR `Athena++ applies source terms at BOTH ...` → DELETE
- L899-901 OVERJUSTIFY → `2. Gas rho/momentum are taken from the post-flux OneCell[]`
- L902-904 OVERJUSTIFY → `3. Applies the gas back-reaction when CR_source is on`
- L905-908 OVERJUSTIFY/DECOR (PRE-stage w ...) → `4. ec_source uses grad(Pc) from the half-step CR fluxes and the t^n gas velocity (Athena++ evaluates it from the pre-stage w)`
- L909-913 OVERJUSTIFY/DECOR (NOT the CT-half-step B ... full-step clause) → `5. Rotation uses the t^n B_n[], not the CT-updated OneCell[MAG_OFFSET+*]` (renumber "see Note 7")
- L938-939 STALE `already done in the main Hydro_RiemannPredict loop above where out_con is updated` → DELETE
- L945 RESTATE `Half-step uses dt_source = 0.5 * dt (matching Athena++ VL2 stage 1 with beta=0.5)` → DELETE
- L948-953 OVERJUSTIFY 6-line CR_E floor story → `1. get the CR state` + `--> floor CR_E here since the passive-scalar floor in Hydro_RiemannPredict() is applied after this function`
- L959-963 OVERJUSTIFY/HISTORY/DECOR (STAGE-UPDATED ... "floored cells will not match") → `2. get gas density and velocity` + `--> floor rho here since the density floor is applied after this function`
- L969-971 OVERJUSTIFY/DECOR (NOT the CT-half-step ... MHD compile-enforced ...) → `3. get the t^n B field (see Note 5)`
- L980-981 STALE/HISTORY/DECOR (READ ... Python/Athena++ pattern) → `5. get sigma_adv and v_adv (set by CR_UpdateStreaming())`
- L991 RESTATE → `total velocity = gas + streaming velocity`
- L1015 STALE → `perpendicular work is added via ec_source (step 11)`
- L1020 RESTATE → `8. effective opacity 1/sigma = 1/sigma_diff + 1/sigma_adv (streaming only)`
- L1035-1036 DECOR `·` → `*`
- L1073-1075 OVERJUSTIFY → `CR_Ec_source=0: drop the Ec source term --> dEc/dt = -div(Fc) (Jiang & Oh 2018 Eq. 19)`
- L1088-1093 OVERJUSTIFY/DECOR (FULL divergence ... 9 terms total) → `11. perpendicular heating ec_source = v_perp*grad(Pc)` + `--> grad(Pc)[n] = div(flux of CR_F(n+1))/vmax, as in CR_UpdateStreaming()`
- L1096 DECOR `CR_F* field indices DESCEND, hence the explicit list` → `CR_F1/2/3 are in descending order`
- L1108-1110 OVERJUSTIFY/DECOR → `perpendicular gas velocity in the B frame (t^n, see Note 4)`
- L1125-1127 OVERJUSTIFY → `13. back-reaction to gas momentum and energy`

### CR_TwoMomentSource_FullStep (near-copy of HalfStep → same fixes)
- L1157-1158 RESTATE `Solves the coupled source terms ... / Includes rotation ...` → `1. Invoked by CPU_FluidSolver_MHM() after Hydro_FullStepUpdate()`
- L1160-1163 OVERJUSTIFY/HISTORY("still") → `3. Gas rho/momentum are taken from the post-flux g_Output[]; B and ADV_* from g_PriVar_Half[]`
- L1164-1168 OVERJUSTIFY/DECOR → `4. ec_source uses grad(Pc) from the full-step CR fluxes and the half-step gas velocity (Athena++ evaluates it from the pre-stage w)`
- L1169 RESTATE `Based on Athena++ cr_source.cpp` → DELETE (Reference says it)
- L1213-1217 OVERJUSTIFY (FMAX "cheap defensive floor so this function stays correct on its own") → `1. get the CR state` + `--> CR_E is already floored by Hydro_FullStepUpdate(); FMAX is only a safeguard`
- L1223-1228 OVERJUSTIFY/HISTORY/DECOR → `2. get gas density and velocity` + `--> floor rho since Hydro_FullStepUpdate() does not (deferred to Flu_Close())`
- L1235 OVERJUSTIFY `(MHD is compile-enforced for CR_TWOMOMENT, see Aux_Check_Parameter.cpp)` → DELETE
- L1244-1245 STALE/HISTORY/DECOR → `5. get sigma_adv and v_adv (set by CR_UpdateStreaming())`
- L1252 RESTATE → `total velocity = gas + streaming velocity`
- L1262 STALE `save original CR flux for back-reaction` → `6. save the original CR state for the floor and back-reaction`
- L1276 STALE → same as L1015
- L1281 RESTATE → same as L1020
- L1296-1299 DECOR `·` → `*`; matrix-form line optional
- L1327-1328 RESTATE → DELETE, keep L1326 (matches half-step L1063)
- L1339-1341 OVERJUSTIFY → same as L1073
- L1354-1360 OVERJUSTIFY/DECOR → same as L1088
- L1364-1366 OVERJUSTIFY → `flux index` + `--> one additional flux is evaluated along each transverse direction for the CT electric field` (upstream CR_AdiabaticWork_FullStep wording)
- L1369 DECOR → same as L1096
- L1381-1385 OVERJUSTIFY/DECOR (`Note: primitive array layout is [DENS=0, VelX=1, ...]`, deliberately ... PRE-stage w) → `perpendicular gas velocity in the B frame (half-step, see Note 4)`

---------------------------------------------------------------------------------------------------
## 2. src/Model_Hydro/CPU_Hydro/CPU_FluidSolver_MHM.cpp
- L42, L159 OVERJUSTIFY/DECOR `CR_TWOMOMENT ... is a standalone module and does NOT require COSMIC_RAY` → DELETE both
- L244-246 PARAM `--> non-const under CR_TWOMOMENT: ... no other ghost-padded array persisting through the full step is available` → `--> Non-const for CR_TWOMOMENT since ADV_* are updated in place`
- L488-491 RESTATE/DECOR (CENTRAL DIFFERENCE, `±1`, output:/input: narration) → `// update ADV_* in place using central-difference grad(Pc)`
- L514-515 RESTATE/HISTORY (`(DefaultStreaming)`, output:/input:) → `// update ADV_* using grad(Pc) from the half-step CR fluxes`
- L564-568 OVERJUSTIFY (retry restore, "regenerates bitwise ...") → `// restore ADV_* overwritten by CR_UpdateStreaming() in the previous iteration`
- L643-644 RESTATE/HISTORY → `// update ADV_* using grad(Pc) from the full-step CR fluxes`
- L670-672 OVERJUSTIFY (widens FLU_GHOST_SIZE by one) → `// > LR_GHOST_SIZE for CR_TWOMOMENT (see FLU_GHOST_SIZE in Macro.h)` or DELETE
- L702-704 HISTORY/DECOR/STALE `NOTE: no opacity update here; Athena++'s end-of-stage-2 DefaultOpacity ... (step 1-a-1)` → DELETE
- L941-942 PARAM → `g_FC_B_In : Array storing the input face-centered B field (for CR_TWOMOMENT only)`
- L1029-1031 OVERJUSTIFY (`deliberately does NOT read out_con[MAG...]`) → DELETE (and the B block can move back to its upstream spot → smaller diff)
- L1044-1051 HISTORY/OVERJUSTIFY (stage-1 b_angle, CRTC_OPACITY, bitwise identical, MHD compile-enforced) →
  `//    t^n cell-centered B for the CR source term`
  `//    --> recomputed since g_PriVar[MAG_OFFSET+*] aliases g_PriVar_Half[], which this loop overwrites`
- L1096-1101 STALE/RESTATE/DECOR (`(DefaultOpacity)`, CENTRAL DIFFERENCE, "next step's flux", NOTE: PRIMITIVE ...) → `// update ADV_* of the half-step solution for the full-step CR fluxes`
- L1106 RESTATE `NSize=N_HF_VAR-2 (interior cells), out_offset=1 (skip boundary), in_offset=1` → DELETE

## 3. Other infrastructure files
- configure.py L707-710 STALE/OVERJUSTIFY → `"Enable two-moment cosmic-ray transport (Jiang & Oh 2018). Must enable <--mhd> and use <--eos=GAMMA>.\n"` (and fix constraint on L706)
- Macro.h L321-328 HISTORY/RESTATE/DECOR (REDUCED, "same convention as Athena++'s u_cr(1:3)") →
  `// two-moment CR fields (Jiang & Oh 2018)`
  `//    CR_E      : CR energy density`
  `//    CR_F1/2/3 : reduced CR flux Fc/CR_VMAX`
  `//    ADV_*     : streaming opacity/velocity work arrays (not passive scalars)`
- Microphysics.h L26-27 STALE (yours) → `CR_sigma : Opacities parallel/perpendicular to the magnetic field = CR_VMAX*CR_SIGMA[_PERP]` / `CR_sigma_perp   (runtime parameters: CR_SIGMA, CR_SIGMA_PERP)`
- Aux_Check_Parameter.cpp L1949 OVERJUSTIFY `(CR transport is field-aligned)` → drop (match CR_DIFFUSION check)
- Aux_Check_Parameter.cpp L1964 OVERJUSTIFY `one extra ghost layer for the outermost ADV_* ...` → DELETE (dup of Macro.h L697)
- Aux_Check_Parameter.cpp L1973-1974 (yours) message never names the parameter → `"OPT__1ST_FLUX_CORR is not supported for CR_TWOMOMENT !!\n"`
- Aux_Check_Parameter.cpp L1984-1985 STALE (yours) "begood" → `"WARNING : CR_CFL (%14.7e) > DT__FLUID (%14.7e) --> CR_CFL has no effect when |v|+c_fast > CR_VMAX !!\n"` (or delete)
- Init_Field.cpp L143 SPECULATIVE (yours) → `// ADV_* work arrays (TODO: allocate only when CR_STREAM is on; see gamer-project/gamer#548)`; put `//` at col 0 like L139-141 (also L152)
- Init_Field.cpp L152 STALE → `// CR energy and flux`
- Main.cpp L373 STALE → `c. cosmic-ray two-moment`
- Microphysics_Init.cpp L41-43 OVERJUSTIFY (yours, ungrammatical) → `// CR_SIGMA[_PERP] = sigma' of Jiang & Oh (2018) with D = 1/(3*sigma'); the solver uses CR_VMAX*sigma'`
- Mis_GetTimeStep.cpp L239 DECOR (yours) `(MIGHT BE REDUNDANT WITH CRITERION ONE??)` → `1.10 CRITERION TEN : cosmic-ray reduced speed of light` — or better: delete the criterion (it IS redundant)
- Mis_GetTimeStep.cpp L241-244 STALE/HISTORY → delete with the criterion
- CPU_dtSolver_HydroCFL.cpp L207-215 HISTORY/OVERJUSTIFY (Athena new_blockdt.cpp, FREEZE fallback, parity) →
  `//    CR signal speed = max( |v|+c_fast, CR_VMAX )`
  `//    --> CR_VMAX is uniform, so max() can be applied after the MaxCFL reduction`
- CPU_Shared_DataReconstruction.cpp L213-215 & L731-733 OVERJUSTIFY → `// > LR_GHOST_SIZE for CR_TWOMOMENT (see FLU_GHOST_SIZE in Macro.h)`
- CPU_Shared_RiemannSolver_HLLE.cpp L678-681 OVERJUSTIFY → DELETE (a generic solver shouldn't document a caller's override; HLLC/HLLD/Roe have none) or `--> overwritten by CR_TwoMomentFlux_*() for CR_TWOMOMENT`
- Clean: perlmutter.config, Field.h, Global.h, Typedef.h, Aux_TakeNote, CPU_FluidSolver, CUAPI_*, Init_Load_Parameter, Init_ResetParameter, Init_TestProb, Makefile_base, CUFLU_CR_TwoMoment.cu

---------------------------------------------------------------------------------------------------
## 4. src/TestProblem/Hydro/CR_Streaming/Init_TestProb_Hydro_CR_Streaming.cpp

### file header / globals
- L5-10 STALE/OVERJUSTIFY/DECOR (NOT set here, "each example directory") → `CR_Streaming_Test options (Jiang & Oh 2018, Sec 4)`
- L18-19 OVERJUSTIFY → `// ports of the classic (COSMIC_RAY) tests CR_Diffusion/ShockTube/SoundWave`
- L20 STALE `frozen gas (pure diffusion)` → `Gaussian Ec, uniform B||x`
- L23, L25 RESTATE section banners → DELETE (L24/L26 say it)
- L32 STALE → `streaming/wave axis (0/1/2 = x/y/z; 3 = (1,1,1), CLASSIC_SOUNDWAVE only)`
- L33-35 OVERJUSTIFY → `background flow speed along CR_Streaming_Dir; inflow speed for CR_TEST_PAPER_SHOCK`
- L37-39 OVERJUSTIFY → `CR_TEST_PAPER_SHOCK initial Fc: 1 = advective (4/3)v*Ec, 0 = zero`
- L40-44 OVERJUSTIFY/RESTATE → `CR_TEST_BLAST_3D: polar angle of B from +z (deg)` / `CR_TEST_BLAST_3D: azimuth of B from +x (deg)`; DELETE L42-44
- L45-49 OVERJUSTIFY/HISTORY ("public-Athena outflow copy") → `1 = gradient-preserving CR outflow BC on the x faces (see GradOutflowBC())`

### Validate / SetParameter
- L98-99 OVERJUSTIFY/FILLER → DELETE (error message suffices)
- L118 STALE template line → DELETE (same at L829)
- L168-169 OVERJUSTIFY `guard against silently ignoring ...` → DELETE
- L212-213 OVERJUSTIFY `... so enabling the flag is a single switch ...` → `CR_Streaming_GradOutflowBC requires the user BC on the x faces`

### CR_Bottleneck_Dens
- L256 DECOR `Uses ABSOLUTE coordinates` → `x is absolute; assumes the paper domain x in [0,1000]`
- L257-258 RESTATE (formula repeats code; Eq. 24 already cited) → DELETE

### SetGridIC
- L287-289 RESTATE/FILLER (`--> we therefore set safe no-streaming defaults`) → `5. ADV_* are overwritten by CR_UpdateOpacity() before use`
- L344-351 OVERJUSTIFY/HISTORY (8-line implosion story, "rho_center 1 -> 2.7 by t = 0.26", "< 1e-5 cells") →
  `Sec 4.1.5: Ec = 12 in 0.5 < r < 0.7 and |phi| < pi/12, else 10; circular B (SetBFieldIC)`
  `rho = 1e8 emulates the paper's frozen MHD variables (circular B has unbalanced hoop stress)`
- L356 RESTATE trailing `quasi-static gas (see the note above)` → DELETE
- L362-363 RESTATE → `Sec 4.1.3: cold cloud; CRs injected from -x by BottleneckBC()`
- L364 RESTATE trailing → DELETE
- L388-391 OVERJUSTIFY → `3D version of Sec 4.2.3: spherical Ec = 100 for r < 0.02, else 0.1; B direction from CR_Streaming_B_theta/B_phi`
- L401-410 OVERJUSTIFY/DECOR/HISTORY/STALE (EXACT classic main-branch setup ...) →
  `Classic CR_Diffusion (Type 0, Mag_Type 0): Ec = E0*exp(-R02*r^2) + BG, uniform B||x, rho = 1, Pgas = 5/3`
  `--> set CR_SIGMA = 1/(3*kappa_classic); CR_SOURCE = 0 keeps the gas static (pure diffusion)`
- L420-423 RESTATE (L/R states repeat code) → `Classic CR_ShockTube (gamma_cr = 4/3); requires large CR_SIGMA, CR_STREAM = 0, CR_SOURCE = 1`
- L434-438 RESTATE → `Classic CR_SoundWave; CR settings as CR_TEST_CLASSIC_SHOCKTUBE; CR_Streaming_Dir = 3 -> along (1,1,1)`
- L465-475 OVERJUSTIFY/DECOR/HISTORY/STALE (paper's OWN ... NOT the classic ... NOTE: REDUCED ...) →
  `Sec 4.2.2 (Fig 11): streams at +/-FlowV collide at xc; rho = Pgas = 1, Ec = CR_Streaming_Ec0, B||x`
  `paper: FlowV = 10, Ec0 = 1/50/200, CR_SIGMA = 10, CR_VMAX = 100`
  `CR_F stores the reduced flux Fc/CR_VMAX`
- L494 FILLER trailing `using the gamma-law EoS, not including the CR pressure` → `gas only`
- L500-501 RESTATE/STALE → `CR fields`
- L507 RESTATE/FILLER → DELETE (repeats Note 5)

### SetBFieldIC
- L532-533 OVERJUSTIFY Note 2 (v_A = |B|/sqrt(rho) ...) → DELETE
- L566-576 HISTORY/DECOR/OVERJUSTIFY (DISCRETE CURL, MACHINE PRECISION, "the previous pointwise field ...", "~3e-8") → `circular B = (-(y-yc), x-xc)/r (Sec 4.1.5), set as the discrete curl of Az = -r so that div(B) = 0 to round-off`
- L579-581 FILLER trailing (`at the correct face center`) → `each component is evaluated at its own face center`
- L591-592 RESTATE → `B || wave vector`
- L606-608 RESTATE → `uniform |B| = 1 along (theta, phi) in degrees; default (90, 0) = B||x as in CR_TEST_BLAST_2D`
- L620-622 RESTATE trailing (×3 "B||x ...") → DELETE (L624 says it)

### BottleneckBC
- L642 FILLER → `2. Requires OPT__BC_FLU_XM = 4`
- L643-650 HISTORY/OVERJUSTIFY ("~2% high relative to Athena++ (the only place the two codes disagreed)") → `3. Ec_ghost = 3, Fc_ghost = -Fc_active (reflecting), following Jiang & Oh 2018 Sec 4.1.3`
- L651-652 OVERJUSTIFY (JetICMWall reference) → `4. Fc_active is read from Array, whose interior is already filled`
- L673 RESTATE → `start from the IC`
- L676 RESTATE → `Ec reservoir`
- L681-683 RESTATE → `reflecting Fc: Fc_ghost = -Fc_active`
- L687 RESTATE trailing → `first active cell`

### ShockBC
- L703-705 OVERJUSTIFY → `User BC for CR_TEST_PAPER_SHOCK: fix the x ghost zones to the IC (Sec 4.2.2)`
- L708-709 OVERJUSTIFY → DELETE (or `2. Valid until the shocks reach the boundaries`)
- L711 PARAM → `See BottleneckBC()` (upstream wording)

### GradOutflowBC
- L731-732 DECOR (BOTH) → `1. Linked to BC_User_Ptr; SetParameter() resets OPT__BC_FLU_XM/XP to user`
- L733-735 OVERJUSTIFY/STALE → `2. Gas fields are reset to the IC: valid only for a uniform gas decoupled from CRs (CR_SOURCE = 0)`
- L736-742 HISTORY/OVERJUSTIFY (LINEARLY EXTRAPOLATED, "public-Athena copy BC halves ... flattens Ec over the last ~12 cells") →
  `3. Ec: linear extrapolation from the last two active cells, floored at CR_EC_TINY; Fc: copied from the last active cell`
  `   --> keeps grad(Pc) at the boundary so CRs stream out at (4/3)v_A*Ec (Jiang & Oh 2018 Fig 3)`
- L743-745 RESTATE → DELETE (repeats SetGridIC Note 5)
- L746-747 RESTATE → DELETE
- L749 PARAM → `See BottleneckBC()`
- L758 FILLER → `start from the IC`

### Init_TestProb_Hydro_CR_Streaming
- L829 STALE template line → DELETE
- L835-840 OVERJUSTIFY (optional) "procedure to enable a problem-specific function" → DELETE (upstream CR_Diffusion drops it)
- L846-847 RESTATE → `spherical-shell refinement (see Flag_CR_Streaming())`
- L854-855 OVERJUSTIFY → `uniform B --> reuse SetBFieldIC() for the B ghost zones`
- L870-872 OVERJUSTIFY → `gradient-preserving CR outflow BC on the x faces; B ghost zones from SetBFieldIC()` (note: assumes uniform B, nothing prevents combining with CR_TEST_CIRCLE_2D)

## 5. Flag_CR_Streaming.cpp
- L8 FILLER (British spelling) → `Flag a spherical shell around the box center`
- L13-14 OVERJUSTIFY → DELETE (or `3. Position-only --> static hierarchy`)
- L40 RESTATE → DELETE (repeats Note 4)

## 6. CR_Diffusion/Init_TestProb_Hydro_CR_Diffusion.cpp (branch lines)
- L237-239 HISTORY/STALE (non-existent CR_Classic_Diffusion dir; explains a removed check) → DELETE
- L240-243 OVERJUSTIFY (hoop stress paragraph, Tseng 2024 Fig 19, a = B^2/(r*rho)) → `// Type 1: circular B has unbalanced hoop stress --> need a large Rho0 to keep the gas quasi-static`

## 7. example/test_problem/Hydro/CR_Streaming/
### Input__Parameter
- L136-137 OVERJUSTIFY → DELETE (L129 marks GAMMA_CR "COSMIC_RAY only") or `# --> P_cr = E_cr/3 (GAMMA_CR unused)`
- L138, L139, L149 STALE `[0/1]` → `include CR source terms (CR-gas coupling) [0]` / `enable CR streaming [1]` / `add CR sound speed to v_diff [0]`
- L140-141 OVERJUSTIFY → `include the CR-energy source term (0 -> pure Eq. 19 of Jiang & Oh 2018) [1]`
- L142 FILLER/STALE → `reduced speed of light for CRs [1e2]`
- L143-144 DECOR/HISTORY (PAPER units, internal/Athena-input sigma) → `inverse CR diffusion coefficient along B: kappa = 1/(3*CR_SIGMA) [1e8]`
- L145 FILLER → `CR_SIGMA perpendicular to B [1e8]`
### Input__TestProb
- L11 STALE → `(Gaussian, uniform B||x)`
- L17 FILLER → `test problem (see above) [0]`
- L18 STALE → add `; 3 = (1,1,1), Test 9 only` and `[0]`; L19 missing `[0.0]`
- L20 FILLER → `gradient-preserving CR outflow BC on the x faces (resets OPT__BC_FLU_XM/XP to 4) [0]`
- L24 RESTATE → `Test 11 only: azimuthal angle of B in degrees (from +x axis) [0.0]`

---------------------------------------------------------------------------------------------------
## 8. Non-comment issues spotted in passing (not asked for; for your triage)
- **BUG (verified)**: TestProb L314/L326 read `amr->BoxCenter[CR_Streaming_Dir]`; ReadPara allows Dir=3 (L153) → out-of-bounds read for CR_TEST_TRIANGULAR_1D / GAUSSIAN_1D with Dir=3. Restrict Dir≤2 for those tests in Validate/SetParameter.
- Criterion Ten (Mis_GetTimeStep L239-248) + dtSolver L211-213 fallback: fully redundant → delete.
- CR_UpdateStreaming/CR_UpdateOpacity: NVar_Out/NVar_In/NVar_B and out_offset/in_offset are equal at every call site → API could be simplified.
- `#ifdef MHD` guards in CR_TwoMoment L1010, 1084, 1094, 1271, 1350, 1361 are dead (MHD is #error-enforced).
- clean.sh L6 removes `CR_Diffusion*` (leftover from the CR_Diffusion example).
- MHM L975 `#ifdef COSMIC_RAY` → `#if ( defined COSMIC_RAY )` with no other change (diff noise); HLLE L686/L693 whitespace-only edits.
- DataReconstruction GAMER_DEBUG check `NIn - 2*NGhost != N_FC_VAR` is now always false by construction.
- CR_TwoMomentFlux_HalfStep: `vmax` is loop-invariant (L689) and could be hoisted.

Please identify all of the non-comment changes in the current staged code, report them to me in short. And verify that if they will not cahnge the code behvior at all (assume the compiler works as it should be).
