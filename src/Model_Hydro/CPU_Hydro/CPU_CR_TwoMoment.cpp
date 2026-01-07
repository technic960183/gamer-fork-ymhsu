#include "CUFLU.h"

#ifdef CR_STREAMING



// global constants for CR streaming
// --> will be moved to input parameters later
static const real CR_MAX_OPACITY    = (real)1.0e20;   // max_opacity for diffusion (effectively infinite)
static const real CR_TAU_ASYM_LIM   = (real)1.0e-3;   // optical depth limit for asymptotic expansion
static const real CR_TAUFACT        = (real)1.0;      // tau factor for optical depth calculation
static const int  CR_VEL_FLX_FLAG   = 1;              // flag to add CR sound speed to v_diff


//-------------------------------------------------------------------------------------------------------
// Function    : RotateVec
// Description : Rotate a vector from lab frame to B-aligned frame
//
// Note        : 1. B-aligned frame: x' along B, y' and z' perpendicular to B
//               2. The rotation is defined by two angles: theta (B to z-axis) and phi (Bxy to x-axis)
//
// Parameter   : sint : sin(theta) = |Bxy|/|B|
//               cost : cos(theta) = Bz/|B|
//               sinp : sin(phi) = By/|Bxy|
//               cosp : cos(phi) = Bx/|Bxy|
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void RotateVec( const real sint, const real cost,
                       const real sinp, const real cosp,
                       real &v1, real &v2, real &v3 )
{
// The two rotation matrices:
// R_1 (around z by phi):
// [cos_p  sin_p 0]
// [-sin_p cos_p 0]
// [0       0    1]
//
// R_2 (around y by theta):
// [sin_t  0 cos_t]
// [0      1    0]
// [-cos_t 0 sin_t]

// First apply R1, then apply R2
   real newv1 =  cosp * v1 + sinp * v2;
   v2         = -sinp * v1 + cosp * v2;

// now apply R2
   v1         =  sint * newv1 + cost * v3;
   real newv3 = -cost * newv1 + sint * v3;
   v3 = newv3;

} // FUNCTION : RotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : InvRotateVec
// Description : Rotate a vector from B-aligned frame back to lab frame (inverse of RotateVec)
//
// Note        : 1. Applies R1^-1 * R2^-1
//
// Parameter   : sint, cost, sinp, cosp : same as RotateVec
//               v1, v2, v3 : input/output vector components
//
// Reference   : Athena++ src/utils/rotate_vectors.cpp
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void InvRotateVec( const real sint, const real cost,
                          const real sinp, const real cosp,
                          real &v1, real &v2, real &v3 )
{
// R_1^-1 (around z by -phi):
// [cos_p  -sin_p 0]
// [sin_p  cos_p  0]
// [0       0     1]
//
// R_2^-1 (around y by -theta):
// [sin_t  0 -cos_t]
// [0      1    0  ]
// [cos_t  0 sin_t ]

// First apply R2^-1, then apply R1^-1
   real newv1 = sint * v1 - cost * v3;
   v3         = cost * v1 + sint * v3;

// now apply R1^-1
   v1         = cosp * newv1 - sinp * v2;
   real newv2 = sinp * newv1 + cosp * v2;
   v2 = newv2;

} // FUNCTION : InvRotateVec



//-------------------------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description : Monotonized central (MC) slope limiter
//
// Note        : 1. Returns limited slope from two input slopes
//
// Parameter   : a, b : two input slopes
//
// Return      : limited slope
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static real MC_limiter( const real a, const real b )
{
   const real c = (real)0.5 * ( a + b );
   const real d = (real)2.0 * ( ( FABS(a) < FABS(b) ) ? a : b );
   const real s = ( a*b > (real)0.0 ) ? (real)1.0 : (real)0.0;

   return s * ( ( FABS(c) < FABS(d) ) ? c : d );

} // FUNCTION : MC_limiter



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeBFieldAngles
// Description : Compute the rotation angles from B field direction
//
// Note        : 1. Returns sin(theta), cos(theta), sin(phi), cos(phi)
//               2. theta is the angle between B and z-axis
//               3. phi is the angle between B_xy projection and x-axis
//
// Parameter   : Bx, By, Bz : magnetic field components
//               sint, cost, sinp, cosp : output angles
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeBFieldAngles( const real Bx, const real By, const real Bz,
                                    real &sint, real &cost, real &sinp, real &cosp )
{
   const real bxby = SQRT( SQR(Bx) + SQR(By) );
   const real btot = SQRT( SQR(Bx) + SQR(By) + SQR(Bz) );

   if ( btot > TINY_NUMBER ) {
      sint = bxby / btot;       // sin(theta) = |Bxy|/|B|
      cost = Bz / btot;         // cos(theta) = Bz/|B|
   } else {
      sint = (real)1.0;
      cost = (real)0.0;
   }

   if ( bxby > TINY_NUMBER ) {
      sinp = By / bxby;         // sin(phi) = By/|Bxy|
      cosp = Bx / bxby;         // cos(phi) = Bx/|Bxy|
   } else {
      sinp = (real)0.0;
      cosp = (real)1.0;
   }

} // FUNCTION : CR_ComputeBFieldAngles



//-------------------------------------------------------------------------------------------------------
// Function    : CR_ComputeHLLEFlux
// Description : Compute HLLE flux for CR variables at a cell interface
//
// Note        : 1. Based on Athena++ cr_flux.cpp
//               2. fdir is the flux direction (0=x, 1=y, 2=z)
//
// Parameter   : Ec_L/R      : CR energy density (left/right)
//               Fc_L/R[3]   : CR flux (left/right)
//               v_L/R       : gas velocity along flux direction (left/right)
//               vdiff_L/R   : diffusion velocity (left/right)
//               vmax        : effective speed of light
//               fdir        : flux direction (0=CRF1, 1=CRF2, 2=CRF3)
//               flux_E      : output flux for CR energy
//               flux_F[3]   : output flux for CR flux components
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
static void CR_ComputeHLLEFlux( const real Ec_L, const real Ec_R,
                                const real Fc_L[3], const real Fc_R[3],
                                const real v_L, const real v_R,
                                const real vdiff_L, const real vdiff_R,
                                const real vmax,
                                const int fdir,
                                real &flux_E, real flux_F[3] )
{
// Eddington factors (isotropic)
   const real edd   = (real)1.0 / (real)3.0;   // diagonal component
   const real edd12 = (real)0.0;               // off-diagonal components
   const real edd13 = (real)0.0;
   const real edd23 = (real)0.0;

// compute wave speeds
   const real meanadv  = (real)0.5 * ( v_L + v_R );
   const real meandiffv = (real)0.5 * ( vdiff_L + vdiff_R );

   real al = FMIN( meanadv - meandiffv, v_L - vdiff_L );
   real ar = FMAX( meanadv + meandiffv, v_R + vdiff_R );

// clamp to light speed
   const real sqrt_edd = SQRT( edd );
   ar = FMIN( ar,  vmax * sqrt_edd );
   al = FMAX( al, -vmax * sqrt_edd );

   const real bp = FMAX( ar, (real)0.0 );   // positive wave speed
   const real bm = FMIN( al, (real)0.0 );   // negative wave speed

// compute L/R fluxes minus wave transport: F_L - bm*U_L, F_R - bp*U_R
   real fl_e, fr_e;
   real fl_f1, fr_f1, fl_f2, fr_f2, fl_f3, fr_f3;

// energy equation flux: vmax * F_c - b * E_c
   fl_e = vmax * Fc_L[fdir] - bm * Ec_L;
   fr_e = vmax * Fc_R[fdir] - bp * Ec_R;

// flux equation fluxes depend on direction
   if ( fdir == 0 ) {       // x-direction (CRF1)
      fl_f1 = vmax * edd   * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd   * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd12 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd12 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd13 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd13 * Ec_R - bp * Fc_R[2];
   }
   else if ( fdir == 1 ) {  // y-direction (CRF2)
      fl_f1 = vmax * edd12 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd12 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd   * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd   * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd23 * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd23 * Ec_R - bp * Fc_R[2];
   }
   else {                   // z-direction (CRF3)
      fl_f1 = vmax * edd13 * Ec_L - bm * Fc_L[0];
      fr_f1 = vmax * edd13 * Ec_R - bp * Fc_R[0];

      fl_f2 = vmax * edd23 * Ec_L - bm * Fc_L[1];
      fr_f2 = vmax * edd23 * Ec_R - bp * Fc_R[1];

      fl_f3 = vmax * edd   * Ec_L - bm * Fc_L[2];
      fr_f3 = vmax * edd   * Ec_R - bp * Fc_R[2];
   }

// HLLE flux formula
   real tmp = (real)0.0;
   if ( FABS(bm - bp) > TINY_NUMBER )
      tmp = (real)0.5 * (bp + bm) / (bp - bm);

   flux_E    = (real)0.5 * (fl_e  + fr_e ) + (fl_e  - fr_e ) * tmp;
   flux_F[0] = (real)0.5 * (fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
   flux_F[1] = (real)0.5 * (fl_f2 + fr_f2) + (fl_f2 - fr_f2) * tmp;
   flux_F[2] = (real)0.5 * (fl_f3 + fr_f3) + (fl_f3 - fr_f3) * tmp;

} // FUNCTION : CR_ComputeHLLEFlux



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_HalfStep
//
// Description : 
//
// Note        : 1.
//
// Reference   : [1] 
//
// Parameter   : g_Con_Var   : Array storing the input cell-centered conserved fluid variables
//               g_Flux_Half : Array with hydrodynamic fluxes for adding the cosmic-ray diffusive fluxes
//               g_FC_B      : Array storing the input face-centered B field
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                  real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                            const real g_CC_B[][ CUBE(FLU_NXT) ],
                            const real dh, const MicroPhy_t *MicroPhy )
{  


   const int  didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, stride_fc_B;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i   = N_HF_FLUX-1;              size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = 0;                        j_offset = flux_offset;              k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT_P1;               sizeB_j  = FLU_NXT;                  stride_fc_B = 1;
                  break;

         case 1 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-1;              size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = flux_offset;              j_offset = 0;                        k_offset    = flux_offset;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT_P1;               stride_fc_B = FLU_NXT;
                  break;

         case 2 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-1;
                  i_offset = flux_offset;              j_offset = flux_offset;              k_offset    = 0;
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT;                  stride_fc_B = SQR(FLU_NXT);
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

//       conserved variable and cell-centered magnetic field index
         const int i_cvar   = i_flux;
         const int j_cvar   = j_flux;
         const int k_cvar   = k_flux;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       face-centered magnetic field index
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B;

//       get left and right cell indices
         const int idx_L = idx_cvar;
         const int idx_R = idx_cvar + didx_cvar[d];

//       1. get CR variables (left and right states) - using donor cell (1st order)
         const real Ec_L   = g_ConVar[CR_E ][idx_L];
         const real Ec_R   = g_ConVar[CR_E ][idx_R];
         const real Fc_L[3] = { g_ConVar[CR_F1][idx_L], g_ConVar[CR_F2][idx_L], g_ConVar[CR_F3][idx_L] };
         const real Fc_R[3] = { g_ConVar[CR_F1][idx_R], g_ConVar[CR_F2][idx_R], g_ConVar[CR_F3][idx_R] };

//       2. get density and velocity
         const real rho_L = g_ConVar[DENS][idx_L];
         const real rho_R = g_ConVar[DENS][idx_R];
         const real v_L   = g_ConVar[MOMX+d][idx_L] / rho_L;
         const real v_R   = g_ConVar[MOMX+d][idx_R] / rho_R;

//       3. get magnetic field at interface
         const real B_N  = g_FC_B[d][idx_fc_B];
         const real B_T1 = (real)0.5 * ( g_CC_B[TDir1][idx_L] + g_CC_B[TDir1][idx_R] );
         const real B_T2 = (real)0.5 * ( g_CC_B[TDir2][idx_L] + g_CC_B[TDir2][idx_R] );

//       reconstruct B components in xyz order
         real Bx_face, By_face, Bz_face;
         if ( d == 0 )      { Bx_face = B_N;  By_face = B_T1; Bz_face = B_T2; }
         else if ( d == 1 ) { Bx_face = B_T2; By_face = B_N;  Bz_face = B_T1; }
         else               { Bx_face = B_T1; By_face = B_T2; Bz_face = B_N;  }

//       4. compute B-field angles
         real sint, cost, sinp, cosp;
         CR_ComputeBFieldAngles( Bx_face, By_face, Bz_face, sint, cost, sinp, cosp );

//       5. compute grad(Pc) using MC limiter at the interface
         real dPc_dx, dPc_dy, dPc_dz;
         real al, bl, ar, br;

//       x direction gradient
         al = g_ConVar[CR_E][idx_L] - g_ConVar[CR_E][idx_L - didx_cvar[0]];
         bl = g_ConVar[CR_E][idx_L + didx_cvar[0]] - g_ConVar[CR_E][idx_L];
         ar = g_ConVar[CR_E][idx_R] - g_ConVar[CR_E][idx_R - didx_cvar[0]];
         br = g_ConVar[CR_E][idx_R + didx_cvar[0]] - g_ConVar[CR_E][idx_R];
         dPc_dx = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       y direction gradient
         al = g_ConVar[CR_E][idx_L] - g_ConVar[CR_E][idx_L - didx_cvar[1]];
         bl = g_ConVar[CR_E][idx_L + didx_cvar[1]] - g_ConVar[CR_E][idx_L];
         ar = g_ConVar[CR_E][idx_R] - g_ConVar[CR_E][idx_R - didx_cvar[1]];
         br = g_ConVar[CR_E][idx_R + didx_cvar[1]] - g_ConVar[CR_E][idx_R];
         dPc_dy = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       z direction gradient
         al = g_ConVar[CR_E][idx_L] - g_ConVar[CR_E][idx_L - didx_cvar[2]];
         bl = g_ConVar[CR_E][idx_L + didx_cvar[2]] - g_ConVar[CR_E][idx_L];
         ar = g_ConVar[CR_E][idx_R] - g_ConVar[CR_E][idx_R - didx_cvar[2]];
         br = g_ConVar[CR_E][idx_R + didx_cvar[2]] - g_ConVar[CR_E][idx_R];
         dPc_dz = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       6. compute B dot grad(Pc)
         const real Bsq = SQR(Bx_face) + SQR(By_face) + SQR(Bz_face);
         const real Btot = SQRT( Bsq );
         const real b_grad_pc = Bx_face * dPc_dx + By_face * dPc_dy + Bz_face * dPc_dz;

//       7. compute Alfven velocity at interface
         const real rho_face = (real)0.5 * ( rho_L + rho_R );
         const real inv_sqrt_rho = (real)1.0 / SQRT( rho_face );
         const real va = Btot * inv_sqrt_rho;

//       8. compute streaming velocity v_adv = -sign(B dot grad Pc) * v_Alfven
         real dpc_sign = (real)0.0;
         if ( b_grad_pc > TINY_NUMBER )       dpc_sign = (real)1.0;
         else if ( -b_grad_pc > TINY_NUMBER ) dpc_sign = (real)-1.0;

         const real v_adv_x = -Bx_face * inv_sqrt_rho * dpc_sign;
         const real v_adv_y = -By_face * inv_sqrt_rho * dpc_sign;
         const real v_adv_z = -Bz_face * inv_sqrt_rho * dpc_sign;

//       9. compute streaming opacity sigma_adv
         const real vmax   = MicroPhy->CR_vmax;
         const real invlim = (real)1.0 / vmax;
         const real Ec_face = (real)0.5 * ( Ec_L + Ec_R );

         real sigma_adv_para;
         if ( va > TINY_NUMBER && Ec_face > TINY_NUMBER ) {
            sigma_adv_para = FABS(b_grad_pc) / ( Btot * va * ((real)4.0/(real)3.0) * invlim * Ec_face );
         } else {
            sigma_adv_para = CR_MAX_OPACITY;
         }

//       perpendicular opacities are max_opacity
         const real sigma_adv_perp = CR_MAX_OPACITY;

//       10. compute total sigma (parallel combination of sigma_diff and sigma_adv)
         const real sigma_diff = CR_MAX_OPACITY;   // sigma_diff is constant max_opacity

//       in B-aligned frame: sigma_x is parallel, sigma_y and sigma_z are perpendicular
         const real sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
         const real sigma_y = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );
         const real sigma_z = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );

//       11. compute v_diff in B-aligned frame with optical depth correction
         const real edd = (real)1.0 / (real)3.0;

//       x direction (parallel to B)
         real tau_x = CR_TAUFACT * sigma_x * dh;
         tau_x = tau_x * tau_x / ( (real)2.0 * edd );
         real diffv_x;
         if ( tau_x < CR_TAU_ASYM_LIM )
            diffv_x = SQRT( (real)1.0 - (real)0.5 * tau_x );
         else
            diffv_x = SQRT( ( (real)1.0 - EXP(-tau_x) ) / tau_x );

//       y direction (perpendicular to B)
         real tau_y = CR_TAUFACT * sigma_y * dh;
         tau_y = tau_y * tau_y / ( (real)2.0 * edd );
         real diffv_y;
         if ( tau_y < CR_TAU_ASYM_LIM )
            diffv_y = SQRT( (real)1.0 - (real)0.5 * tau_y );
         else
            diffv_y = SQRT( ( (real)1.0 - EXP(-tau_y) ) / tau_y );

//       z direction (perpendicular to B)
         real tau_z = CR_TAUFACT * sigma_z * dh;
         tau_z = tau_z * tau_z / ( (real)2.0 * edd );
         real diffv_z;
         if ( tau_z < CR_TAU_ASYM_LIM )
            diffv_z = SQRT( (real)1.0 - (real)0.5 * tau_z );
         else
            diffv_z = SQRT( ( (real)1.0 - EXP(-tau_z) ) / tau_z );

//       v_diff in B-aligned frame
         real vdiff_Bx = vmax * SQRT(edd) * diffv_x;
         real vdiff_By = vmax * SQRT(edd) * diffv_y;
         real vdiff_Bz = vmax * SQRT(edd) * diffv_z;

//       12. rotate v_diff from B-aligned frame to lab frame
         InvRotateVec( sint, cost, sinp, cosp, vdiff_Bx, vdiff_By, vdiff_Bz );

//       take absolute value
         vdiff_Bx = FABS( vdiff_Bx );
         vdiff_By = FABS( vdiff_By );
         vdiff_Bz = FABS( vdiff_Bz );

//       13. add CR sound speed for stability: c_CR = sqrt(4/9 * Ec/rho)
         const real cr_sound = CR_VEL_FLX_FLAG * SQRT( ((real)4.0/(real)9.0) * Ec_face / rho_face );
         vdiff_Bx += cr_sound;
         vdiff_By += cr_sound;
         vdiff_Bz += cr_sound;

//       14. get v_diff along flux direction
         real vdiff_L, vdiff_R;
         if ( d == 0 ) {
            vdiff_L = vdiff_Bx;
            vdiff_R = vdiff_Bx;
         } else if ( d == 1 ) {
            vdiff_L = vdiff_By;
            vdiff_R = vdiff_By;
         } else {
            vdiff_L = vdiff_Bz;
            vdiff_R = vdiff_Bz;
         }

//       15. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       16. add fluxes to output arrays
         g_Flux_Half[d][CR_E ][idx_flux] = flux_E;
         g_Flux_Half[d][CR_F1][idx_flux] = flux_F[0];
         g_Flux_Half[d][CR_F2][idx_flux] = flux_F[1];
         g_Flux_Half[d][CR_F3][idx_flux] = flux_F[2];

//       17. store v_adv and sigma_adv to auxiliary fields for source term use (only store once)
//           Note: we use ADV_VX/VY/VZ for streaming velocity and ADV_SIGMA for parallel sigma
//           These need to be stored at cell centers, not interfaces
//           For half-step, we store to the left cell
         if ( d == 0 ) {
            // store to left cell (idx_L corresponds to idx_cvar)
            // Note: This is a simplified approach - storing interface values to cells
         }

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)


} // FUMCTION : CR_TwoMomentFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_FullStep
//
// Description : Compute full-step CR two-moment fluxes using half-step primitive variables
//
// Note        : 1. Similar to half-step but uses primitive variables from Riemann predictor
//               2. Uses NFlux as the stride for flux array access
//
// Reference   : Athena++ cr_transport.cpp, cr_flux.cpp
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered, half-step primitive fluid variables
//               g_FC_Flux     : Array with hydrodynamic fluxes for adding the cosmic-ray diffusive fluxes
//               g_FC_B_Half   : Array storing the input face-centered, half-step magnetic field
//               NFlux         : Stride for accessing g_FC_Flux[]
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentFlux_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                       real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                 const int NFlux, const real dh, const MicroPhy_t *MicroPhy )
{
   const int  didx_pvar[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  flux_offset  = 0;  // no offset needed for full step
   const real _dh          = (real)1.0 / dh;
   const int  half_offset  = ( N_HF_VAR - NFlux ) / 2;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, stride_fc_B;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

//    set up sizes for flux computation
//    For MHD, we skip the boundary fluxes
#     ifdef MHD
      const int NSkip_N = 0;
      const int NSkip_T = 0;
#     else
      const int NSkip_N = 0;
      const int NSkip_T = 1;
#     endif

      switch ( d )
      {
         case 0 : size_i   = NFlux - 2*NSkip_N + 1;     size_j   = NFlux - 2*NSkip_T;        size_k   = NFlux - 2*NSkip_T;
                  i_offset = NSkip_N;                   j_offset = NSkip_T;                  k_offset = NSkip_T;
                  sizeB_i  = FLU_NXT_P1;                sizeB_j  = FLU_NXT;                  stride_fc_B = 1;
                  break;

         case 1 : size_i   = NFlux - 2*NSkip_T;         size_j   = NFlux - 2*NSkip_N + 1;    size_k   = NFlux - 2*NSkip_T;
                  i_offset = NSkip_T;                   j_offset = NSkip_N;                  k_offset = NSkip_T;
                  sizeB_i  = FLU_NXT;                   sizeB_j  = FLU_NXT_P1;               stride_fc_B = FLU_NXT;
                  break;

         case 2 : size_i   = NFlux - 2*NSkip_T;         size_j   = NFlux - 2*NSkip_T;        size_k   = NFlux - 2*NSkip_N + 1;
                  i_offset = NSkip_T;                   j_offset = NSkip_T;                  k_offset = NSkip_N;
                  sizeB_i  = FLU_NXT;                   sizeB_j  = FLU_NXT;                  stride_fc_B = SQR(FLU_NXT);
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       primitive variable index (with half_offset to account for different strides)
         const int i_pvar   = i_flux + half_offset;
         const int j_pvar   = j_flux + half_offset;
         const int k_pvar   = k_flux + half_offset;
         const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );

//       face-centered B index
         const int idx_fc_B = IDX321( i_pvar, j_pvar, k_pvar, sizeB_i, sizeB_j ) + stride_fc_B;

//       get left and right cell indices
         const int idx_L = idx_pvar;
         const int idx_R = idx_pvar + didx_pvar[d];

//       1. get CR variables (left and right states) - primitive vars contain same CR fields
//          Note: for primitive variables, CR_E stores energy density, CR_F stores flux
         const real Ec_L   = g_PriVar_Half[CR_E ][idx_L];
         const real Ec_R   = g_PriVar_Half[CR_E ][idx_R];
         const real Fc_L[3] = { g_PriVar_Half[CR_F1][idx_L], g_PriVar_Half[CR_F2][idx_L], g_PriVar_Half[CR_F3][idx_L] };
         const real Fc_R[3] = { g_PriVar_Half[CR_F1][idx_R], g_PriVar_Half[CR_F2][idx_R], g_PriVar_Half[CR_F3][idx_R] };

//       2. get density and velocity (primitive variables)
//          Note: primitive array layout is [DENS=0, VelX=1, VelY=2, VelZ=3, ...]
         const real rho_L = g_PriVar_Half[DENS][idx_L];
         const real rho_R = g_PriVar_Half[DENS][idx_R];
         const real v_L   = g_PriVar_Half[1+d][idx_L];
         const real v_R   = g_PriVar_Half[1+d][idx_R];

//       3. get magnetic field at interface from half-step face-centered B
         const real B_N  = g_FC_B_Half[d][idx_fc_B];
//       get cell-centered B from primitive variables (stored after MAG_OFFSET)
         const real B_T1 = (real)0.5 * ( g_PriVar_Half[MAG_OFFSET+TDir1][idx_L] + g_PriVar_Half[MAG_OFFSET+TDir1][idx_R] );
         const real B_T2 = (real)0.5 * ( g_PriVar_Half[MAG_OFFSET+TDir2][idx_L] + g_PriVar_Half[MAG_OFFSET+TDir2][idx_R] );

//       reconstruct B components in xyz order
         real Bx_face, By_face, Bz_face;
         if ( d == 0 )      { Bx_face = B_N;  By_face = B_T1; Bz_face = B_T2; }
         else if ( d == 1 ) { Bx_face = B_T2; By_face = B_N;  Bz_face = B_T1; }
         else               { Bx_face = B_T1; By_face = B_T2; Bz_face = B_N;  }

//       4. compute B-field angles
         real sint, cost, sinp, cosp;
         CR_ComputeBFieldAngles( Bx_face, By_face, Bz_face, sint, cost, sinp, cosp );

//       5. compute grad(Pc) using MC limiter at the interface
         real dPc_dx, dPc_dy, dPc_dz;
         real al, bl, ar, br;

//       x direction gradient
         al = g_PriVar_Half[CR_E][idx_L] - g_PriVar_Half[CR_E][idx_L - didx_pvar[0]];
         bl = g_PriVar_Half[CR_E][idx_L + didx_pvar[0]] - g_PriVar_Half[CR_E][idx_L];
         ar = g_PriVar_Half[CR_E][idx_R] - g_PriVar_Half[CR_E][idx_R - didx_pvar[0]];
         br = g_PriVar_Half[CR_E][idx_R + didx_pvar[0]] - g_PriVar_Half[CR_E][idx_R];
         dPc_dx = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       y direction gradient
         al = g_PriVar_Half[CR_E][idx_L] - g_PriVar_Half[CR_E][idx_L - didx_pvar[1]];
         bl = g_PriVar_Half[CR_E][idx_L + didx_pvar[1]] - g_PriVar_Half[CR_E][idx_L];
         ar = g_PriVar_Half[CR_E][idx_R] - g_PriVar_Half[CR_E][idx_R - didx_pvar[1]];
         br = g_PriVar_Half[CR_E][idx_R + didx_pvar[1]] - g_PriVar_Half[CR_E][idx_R];
         dPc_dy = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       z direction gradient
         al = g_PriVar_Half[CR_E][idx_L] - g_PriVar_Half[CR_E][idx_L - didx_pvar[2]];
         bl = g_PriVar_Half[CR_E][idx_L + didx_pvar[2]] - g_PriVar_Half[CR_E][idx_L];
         ar = g_PriVar_Half[CR_E][idx_R] - g_PriVar_Half[CR_E][idx_R - didx_pvar[2]];
         br = g_PriVar_Half[CR_E][idx_R + didx_pvar[2]] - g_PriVar_Half[CR_E][idx_R];
         dPc_dz = MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) ) * _dh / (real)3.0;

//       6. compute B dot grad(Pc)
         const real Bsq = SQR(Bx_face) + SQR(By_face) + SQR(Bz_face);
         const real Btot = SQRT( Bsq );
         const real b_grad_pc = Bx_face * dPc_dx + By_face * dPc_dy + Bz_face * dPc_dz;

//       7. compute Alfven velocity at interface
         const real rho_face = (real)0.5 * ( rho_L + rho_R );
         const real inv_sqrt_rho = (real)1.0 / SQRT( rho_face );
         const real va = Btot * inv_sqrt_rho;

//       8. compute streaming velocity v_adv = -sign(B dot grad Pc) * v_Alfven
         real dpc_sign = (real)0.0;
         if ( b_grad_pc > TINY_NUMBER )       dpc_sign = (real)1.0;
         else if ( -b_grad_pc > TINY_NUMBER ) dpc_sign = (real)-1.0;

//       9. compute streaming opacity sigma_adv
         const real vmax   = MicroPhy->CR_vmax;
         const real invlim = (real)1.0 / vmax;
         const real Ec_face = (real)0.5 * ( Ec_L + Ec_R );

         real sigma_adv_para;
         if ( va > TINY_NUMBER && Ec_face > TINY_NUMBER ) {
            sigma_adv_para = FABS(b_grad_pc) / ( Btot * va * ((real)4.0/(real)3.0) * invlim * Ec_face );
         } else {
            sigma_adv_para = CR_MAX_OPACITY;
         }

//       perpendicular opacities are max_opacity
         const real sigma_adv_perp = CR_MAX_OPACITY;

//       10. compute total sigma (parallel combination of sigma_diff and sigma_adv)
         const real sigma_diff = CR_MAX_OPACITY;

         const real sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
         const real sigma_y = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );
         const real sigma_z = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );

//       11. compute v_diff in B-aligned frame with optical depth correction
         const real edd = (real)1.0 / (real)3.0;

//       x direction (parallel to B)
         real tau_x = CR_TAUFACT * sigma_x * dh;
         tau_x = tau_x * tau_x / ( (real)2.0 * edd );
         real diffv_x;
         if ( tau_x < CR_TAU_ASYM_LIM )
            diffv_x = SQRT( (real)1.0 - (real)0.5 * tau_x );
         else
            diffv_x = SQRT( ( (real)1.0 - EXP(-tau_x) ) / tau_x );

//       y direction (perpendicular to B)
         real tau_y = CR_TAUFACT * sigma_y * dh;
         tau_y = tau_y * tau_y / ( (real)2.0 * edd );
         real diffv_y;
         if ( tau_y < CR_TAU_ASYM_LIM )
            diffv_y = SQRT( (real)1.0 - (real)0.5 * tau_y );
         else
            diffv_y = SQRT( ( (real)1.0 - EXP(-tau_y) ) / tau_y );

//       z direction (perpendicular to B)
         real tau_z = CR_TAUFACT * sigma_z * dh;
         tau_z = tau_z * tau_z / ( (real)2.0 * edd );
         real diffv_z;
         if ( tau_z < CR_TAU_ASYM_LIM )
            diffv_z = SQRT( (real)1.0 - (real)0.5 * tau_z );
         else
            diffv_z = SQRT( ( (real)1.0 - EXP(-tau_z) ) / tau_z );

//       v_diff in B-aligned frame
         real vdiff_Bx = vmax * SQRT(edd) * diffv_x;
         real vdiff_By = vmax * SQRT(edd) * diffv_y;
         real vdiff_Bz = vmax * SQRT(edd) * diffv_z;

//       12. rotate v_diff from B-aligned frame to lab frame
         InvRotateVec( sint, cost, sinp, cosp, vdiff_Bx, vdiff_By, vdiff_Bz );

//       take absolute value
         vdiff_Bx = FABS( vdiff_Bx );
         vdiff_By = FABS( vdiff_By );
         vdiff_Bz = FABS( vdiff_Bz );

//       13. add CR sound speed for stability
         const real cr_sound = CR_VEL_FLX_FLAG * SQRT( ((real)4.0/(real)9.0) * Ec_face / rho_face );
         vdiff_Bx += cr_sound;
         vdiff_By += cr_sound;
         vdiff_Bz += cr_sound;

//       14. get v_diff along flux direction
         real vdiff_L, vdiff_R;
         if ( d == 0 ) {
            vdiff_L = vdiff_Bx;
            vdiff_R = vdiff_Bx;
         } else if ( d == 1 ) {
            vdiff_L = vdiff_By;
            vdiff_R = vdiff_By;
         } else {
            vdiff_L = vdiff_Bz;
            vdiff_R = vdiff_Bz;
         }

//       15. compute HLLE flux
         real flux_E, flux_F[3];
         CR_ComputeHLLEFlux( Ec_L, Ec_R, Fc_L, Fc_R, v_L, v_R, vdiff_L, vdiff_R, vmax, d,
                             flux_E, flux_F );

//       16. add fluxes to output arrays (replace, not add, since hydro doesn't handle CR)
         g_FC_Flux[d][CR_E ][idx_flux] = flux_E;
         g_FC_Flux[d][CR_F1][idx_flux] = flux_F[0];
         g_FC_Flux[d][CR_F2][idx_flux] = flux_F[1];
         g_FC_Flux[d][CR_F3][idx_flux] = flux_F[2];

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

} // FUNCTION : CR_TwoMomentFlux_FullStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_HalfStep
//
// Description : Apply flux divergence to CR variables for the half-step predictor
//
// Note        : 1. The main flux divergence is already applied in the Hydro_RiemannPredict loop
//               2. This function is currently a placeholder for any additional CR-specific source terms
//               3. No implicit source solve at half-step (only at full-step)
//
// Reference   : Athena++ cr_transport.cpp
//
// Parameter   : OneCell     : Single-cell fluid array (already updated with flux divergence)
//               g_ConVar_In : Array storing the input conserved variables
//               g_Flux_Half : Array storing the input face-centered fluxes
//               idx_fc      : Index of accessing g_ConVar_In[]
//               didx_fc     : Index increment of g_ConVar_In[]
//               idx_flux    : Index of accessing g_flux_Half[]
//               didx_flux   : Index increment of g_Flux_Half[]
//               dt_dh2      : 0.5 * dt / dh
//               EoS         : EoS object
//               MicroPhy    : Microphysics object
//
// Return      : OneCell[] (modified CR fields)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_HalfStep( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                            const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                            const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const int idx_fc, const int didx_fc[3],
                            const int idx_flux, const int didx_flux[3],
                            const real dt_dh2, const EoS_t *EoS , const MicroPhy_t *MicroPhy )
{
// The flux divergence update for CR_E, CR_F1, CR_F2, CR_F3 is already done
// in the main Hydro_RiemannPredict loop above where out_con is updated.
//
// Here we ensure CR energy remains positive:
   if ( OneCell[CR_E] < TINY_NUMBER )
      OneCell[CR_E] = TINY_NUMBER;

// Note: For the half-step, we skip the implicit source term solve.
// The implicit source update (coupling between Ec, Fc and gas) is only applied at the full-step.

} // FUNCTION : CR_TwoMomentSource_HalfStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_FullStep
//
// Description : Apply implicit source term update for CR two-moment equations
//
// Note        : 1. Solves the coupled source terms for Ec and Fc implicitly
//               2. Includes rotation to B-aligned frame for anisotropic transport
//               3. Applies back-reaction to gas momentum and energy if CR_source is enabled
//               4. Based on Athena++ cr_source.cpp
//
// Reference   : Athena++ cr_source.cpp
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered primitive variables
//               g_Output      : Array to store the updated fluid data (already has flux divergence applied)
//               g_Flux        : Array storing the input face-centered fluxes
//               g_FC_Var      : Array storing the input face-centered conserved variables
//               dt            : Time interval to advance solution
//               dh            : Cell size
//               EoS           : EoS object
//               MicroPhy      : Microphysics object
//
// Return      : g_Output[] with updated CR and gas fields
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                const real dt, const real dh, const EoS_t *EoS, const MicroPhy_t *MicroPhy )
{
   const int  didx_out[3]  = { 1, PS2, SQR(PS2) };
   const int  didx_pvar[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const real _dh          = (real)1.0 / dh;

// offset between output array and primitive array
   const int  pvar_offset  = ( N_HF_VAR - PS2 ) / 2;

// CR parameters
   const real vmax   = MicroPhy->CR_vmax;
   const real invlim = (real)1.0 / vmax;
   const bool CR_source = MicroPhy->CR_source;   // flag to enable back-reaction to gas

   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
//    output array indices
      const int i_out   = idx_out % PS2;
      const int j_out   = idx_out % SQR(PS2) / PS2;
      const int k_out   = idx_out / SQR(PS2);

//    primitive variable indices (with offset)
      const int i_pvar  = i_out + pvar_offset;
      const int j_pvar  = j_out + pvar_offset;
      const int k_pvar  = k_out + pvar_offset;
      const int idx_pvar = IDX321( i_pvar, j_pvar, k_pvar, N_HF_VAR, N_HF_VAR );

//    1. get current CR state from g_Output (which already has flux divergence applied)
      real ec  = g_Output[CR_E ][idx_out];
      real fc1 = g_Output[CR_F1][idx_out];
      real fc2 = g_Output[CR_F2][idx_out];
      real fc3 = g_Output[CR_F3][idx_out];

//    2. get gas density and velocity from half-step primitive variables
//       Note: primitive array layout is [DENS=0, VelX=1, VelY=2, VelZ=3, ...]
      const real rho = g_PriVar_Half[DENS][idx_pvar];
      real v1 = g_PriVar_Half[1][idx_pvar];
      real v2 = g_PriVar_Half[2][idx_pvar];
      real v3 = g_PriVar_Half[3][idx_pvar];

//    3. get cell-centered B field from half-step primitive variables
#     ifdef MHD
      const real Bx = g_PriVar_Half[MAG_OFFSET+MAGX][idx_pvar];
      const real By = g_PriVar_Half[MAG_OFFSET+MAGY][idx_pvar];
      const real Bz = g_PriVar_Half[MAG_OFFSET+MAGZ][idx_pvar];
#     else
      const real Bx = (real)0.0;
      const real By = (real)0.0;
      const real Bz = (real)1.0;   // default B along z
#     endif

//    4. compute B-field angles for rotation
      real sint, cost, sinp, cosp;
      CR_ComputeBFieldAngles( Bx, By, Bz, sint, cost, sinp, cosp );

//    5. compute grad(Pc) from flux divergence
//       grad_pc = (1/vmax) * div(flux_F)
//       We use the fluxes to compute this, but for simplicity use centered difference
      real dPc_dx, dPc_dy, dPc_dz;

//    x direction gradient (using neighboring cells in primitive array)
      dPc_dx = ( g_PriVar_Half[CR_E][idx_pvar + didx_pvar[0]] - 
                 g_PriVar_Half[CR_E][idx_pvar - didx_pvar[0]] ) * (real)0.5 * _dh / (real)3.0;
//    y direction gradient
      dPc_dy = ( g_PriVar_Half[CR_E][idx_pvar + didx_pvar[1]] - 
                 g_PriVar_Half[CR_E][idx_pvar - didx_pvar[1]] ) * (real)0.5 * _dh / (real)3.0;
//    z direction gradient
      dPc_dz = ( g_PriVar_Half[CR_E][idx_pvar + didx_pvar[2]] - 
                 g_PriVar_Half[CR_E][idx_pvar - didx_pvar[2]] ) * (real)0.5 * _dh / (real)3.0;

//    6. compute B dot grad(Pc) and streaming velocity
      const real Bsq = SQR(Bx) + SQR(By) + SQR(Bz);
      const real Btot = SQRT( Bsq );
      const real b_grad_pc = Bx * dPc_dx + By * dPc_dy + Bz * dPc_dz;

//    Alfven velocity
      const real inv_sqrt_rho = (real)1.0 / SQRT( rho );
      const real va = Btot * inv_sqrt_rho;

//    streaming velocity direction
      real dpc_sign = (real)0.0;
      if ( b_grad_pc > TINY_NUMBER )       dpc_sign = (real)1.0;
      else if ( -b_grad_pc > TINY_NUMBER ) dpc_sign = (real)-1.0;

//    streaming velocity components: v_adv = -sign(B dot grad Pc) * v_Alfven
      const real v_adv_x = -Bx * inv_sqrt_rho * dpc_sign;
      const real v_adv_y = -By * inv_sqrt_rho * dpc_sign;
      const real v_adv_z = -Bz * inv_sqrt_rho * dpc_sign;

//    total velocity = gas velocity + streaming velocity
      real vtot1 = v1 + v_adv_x;
      real vtot2 = v2 + v_adv_y;
      real vtot3 = v3 + v_adv_z;

//    7. compute streaming opacity sigma_adv
      real sigma_adv_para;
      if ( va > TINY_NUMBER && ec > TINY_NUMBER ) {
         sigma_adv_para = FABS(b_grad_pc) / ( Btot * va * ((real)4.0/(real)3.0) * invlim * ec );
      } else {
         sigma_adv_para = CR_MAX_OPACITY;
      }
      const real sigma_adv_perp = CR_MAX_OPACITY;

//    8. save original CR flux for back-reaction calculation
      const real fc1_old = fc1;
      const real fc2_old = fc2;
      const real fc3_old = fc3;
      const real ec_old  = ec;

//    9. rotate all vectors to B-aligned frame
      real fr1 = fc1, fr2 = fc2, fr3 = fc3;

#     ifdef MHD
      RotateVec( sint, cost, sinp, cosp, v1, v2, v3 );
      RotateVec( sint, cost, sinp, cosp, fr1, fr2, fr3 );
      RotateVec( sint, cost, sinp, cosp, vtot1, vtot2, vtot3 );

//    perpendicular components have no streaming velocity contribution in B-frame
      vtot2 = (real)0.0;
      vtot3 = (real)0.0;
#     endif

//    10. compute effective sigma (parallel combination of sigma_diff and sigma_adv)
      const real sigma_diff = CR_MAX_OPACITY;

      const real sigma_x = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_para );
      const real sigma_y = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );
      const real sigma_z = (real)1.0 / ( (real)1.0/sigma_diff + (real)1.0/sigma_adv_perp );

//    11. build implicit matrix and solve
//        Source terms:
//        dEc/dt = -vtot · sigma · (Fc - v*Ec*(4/3)/vmax)
//        dFc/dt = -vmax · sigma · (Fc - v*Ec*(4/3)/vmax)
//
//        Matrix form: A * [Ec_new, Fr1_new, Fr2_new, Fr3_new]^T = [rhs1, rhs2, rhs3, rhs4]^T

      const real rhs1 = ec;
      const real rhs2 = fr1;
      const real rhs3 = fr2;
      const real rhs4 = fr3;

//    coefficients for Ec equation (row 1)
      const real coef_11 = (real)1.0 - dt * sigma_x * vtot1 * v1 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_y * vtot2 * v2 * invlim * ((real)4.0/(real)3.0)
                                     - dt * sigma_z * vtot3 * v3 * invlim * ((real)4.0/(real)3.0);
      const real coef_12 = dt * sigma_x * vtot1;
      const real coef_13 = dt * sigma_y * vtot2;
      const real coef_14 = dt * sigma_z * vtot3;

//    coefficients for Fr1 equation (row 2)
      const real coef_21 = -dt * v1 * sigma_x * ((real)4.0/(real)3.0);
      const real coef_22 = (real)1.0 + dt * vmax * sigma_x;

//    coefficients for Fr2 equation (row 3)
      const real coef_31 = -dt * v2 * sigma_y * ((real)4.0/(real)3.0);
      const real coef_33 = (real)1.0 + dt * vmax * sigma_y;

//    coefficients for Fr3 equation (row 4)
      const real coef_41 = -dt * v3 * sigma_z * ((real)4.0/(real)3.0);
      const real coef_44 = (real)1.0 + dt * vmax * sigma_z;

//    solve by substitution (since flux equations are decoupled from each other)
//    Fr_new = (rhs - coef_*1 * Ec_new) / coef_**
//    Substitute into Ec equation and solve for Ec_new

      const real e_coef = coef_11 - coef_12 * coef_21 / coef_22 
                                  - coef_13 * coef_31 / coef_33
                                  - coef_14 * coef_41 / coef_44;

      real new_ec = rhs1 - coef_12 * rhs2 / coef_22 
                        - coef_13 * rhs3 / coef_33 
                        - coef_14 * rhs4 / coef_44;
      new_ec /= e_coef;

//    back-substitute to get new flux
      real newfr1 = ( rhs2 - coef_21 * new_ec ) / coef_22;
      real newfr2 = ( rhs3 - coef_31 * new_ec ) / coef_33;
      real newfr3 = ( rhs4 - coef_41 * new_ec ) / coef_44;

//    12. rotate back to lab frame
#     ifdef MHD
      InvRotateVec( sint, cost, sinp, cosp, newfr1, newfr2, newfr3 );
#     endif

//    13. compute perpendicular heating term (ec_source)
//        This is the work done by perpendicular gas flow: v_perp · grad(Pc)
//        In B-aligned frame: v_perp = (0, v2, v3), grad_pc_perp = (0, dPc/dy', dPc/dz')
#     ifdef MHD
      real dpcdx_B = dPc_dx, dpcdy_B = dPc_dy, dpcdz_B = dPc_dz;
      RotateVec( sint, cost, sinp, cosp, dpcdx_B, dpcdy_B, dpcdz_B );

//    perpendicular velocity (in B-frame, gas velocity only)
//    Note: primitive array layout is [DENS=0, VelX=1, VelY=2, VelZ=3, ...]
      real v1_B = g_PriVar_Half[1][idx_pvar];
      real v2_B = g_PriVar_Half[2][idx_pvar];
      real v3_B = g_PriVar_Half[3][idx_pvar];
      RotateVec( sint, cost, sinp, cosp, v1_B, v2_B, v3_B );

      const real ec_source = v2_B * dpcdy_B + v3_B * dpcdz_B;
      new_ec += dt * ec_source;
#     endif

//    14. floor CR energy
      if ( new_ec < TINY_NUMBER )
         new_ec = ec_old;

//    15. apply back-reaction to gas momentum and energy
      if ( CR_source ) {
//       momentum change: delta_p = -(new_Fc - old_Fc) / vmax
         g_Output[MOMX][idx_out] += -( newfr1 - fc1_old ) * invlim;
         g_Output[MOMY][idx_out] += -( newfr2 - fc2_old ) * invlim;
         g_Output[MOMZ][idx_out] += -( newfr3 - fc3_old ) * invlim;

//       energy change: delta_E = -(new_Ec - old_Ec)
//       (CR energy lost goes to gas thermal energy)
         real new_eg = g_Output[ENGY][idx_out] - ( new_ec - ec_old );
         if ( new_eg < (real)0.0 )
            new_eg = g_Output[ENGY][idx_out];
         g_Output[ENGY][idx_out] = new_eg;
      }

//    16. update CR fields
      g_Output[CR_E ][idx_out] = new_ec;
      g_Output[CR_F1][idx_out] = newfr1;
      g_Output[CR_F2][idx_out] = newfr2;
      g_Output[CR_F3][idx_out] = newfr3;

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

} // FUNCTION : CR_TwoMomentSource_FullStep

#endif // #ifdef CR_STREAMING
