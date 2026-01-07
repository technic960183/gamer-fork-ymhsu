#include "CUFLU.h"

#ifdef CR_STREAMING




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

         //g_Flux_Half[d][CR_E][ idx_flux ] = (real)0.5 * (fl_e + fr_e) + (fl_e - fr_e) * tmp;
         //g_Flux_Half[d][CR_F1][ idx_flux ] = (real)0.5 * (fl_f1 + fr_f1) + (fl_f1 - fr_f1) * tmp;
         
      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

   

   
   // update opacity

} // FUMCTION : CR_TwoMomentFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : CR_TwoMomentFlux_FullStep
//
// Description : 
//
// Note        : 1.
//
// Reference   : [1] 
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

} // FUNCTION : CR_TwoMomentFlux_FullStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_HalfStep
//
// Description : 
//
// Note        : 1.
//
// Reference   : [1] 
//
// Parameter   : OneCell     : Single-cell fluid array to store the updated cell-centered cosmic-ray energy
//               g_ConVar_In : Array storing the input conserved variables
//               g_Flux_Half : Array storing the input face-centered fluxes
//                             --> Accessed with the stride didx_flux
//               idx_in      : Index of accessing g_ConVar_In[]
//               didx_in     : Index increment of g_ConVar_In[]
//               idx_flux    : Index of accessing g_flux_Half[]
//               didx_flux   : Index increment of g_Flux_Half[]
//               dt_dh2      : 0.5 * dt / dh
//               EoS         : EoS object
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_HalfStep( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                            const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                            const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const int idx_fc, const int didx_fc[3],
                            const int idx_flux, const int didx_flux[3],
                            const real dt_dh2, const EoS_t *EoS , const MicroPhy_t *MicroPhy )
{

} // FUMCTION : CR_TwoMomentSource_HalfStep



//-------------------------------------------------------------------------------------------------------
// Function    : CR_TwoMomentSource_FullStep
//
// Description : 
//
// Note        : 1.
//
// Reference   : [1] 
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered primitive variables
//                               --> Accessed with the stride N_HF_VAR
//                               --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//               g_Output      : Array to store the updated fluid data
//               g_Flux        : Array storing the input face-centered fluxes
//                               --> Accessed with the array stride N_FL_FLUX even thought its actually
//                                   allocated size is N_FC_FLUX^3
//               g_FC_Var      : Array storing the input face-centered conserved variables
//                               --> Accessed with the array stride N_FC_VAR^3
//               dt            : Time interval to advance solution
//               dh            : Cell size
//               EoS           : EoS object
//               MicroPhy      : Microphysics object
//
// Return      : g_Output[CRAY][]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_TwoMomentSource_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                const real dt, const real dh, const EoS_t *EoS, const MicroPhy_t *MicroPhy )
{

} // FUMCTION : CR_TwoMomentSource_FullStep

#endif // #ifdef CR_STREAMING
