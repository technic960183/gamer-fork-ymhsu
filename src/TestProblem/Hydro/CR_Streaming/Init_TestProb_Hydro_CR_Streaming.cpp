#include "GAMER.h"



// =======================================================================================
// CR_Streaming_Test : initial-condition geometry selector (Jiang & Oh 2018, Section 4)
// --> the CR *physics* parameters (CR_SIGMA, CR_SIGMA_PERP, CR_STREAM, CR_SOURCE, CR_VMAX, ...)
//     are taken from Input__Parameter, NOT set here.  Each example directory carries the
//     Input__Parameter/Input__TestProb appropriate for one paper test.
// =======================================================================================
#define CR_TEST_TRIANGULAR_1D     0     // Sec 4.1.1 : 1D triangular Ec profile      (streaming)
#define CR_TEST_GAUSSIAN_1D       1     // Sec 4.1.1 / 4.1.4 : 1D Gaussian Ec         (streaming or diffusion)
#define CR_TEST_GAUSSIAN_2D       2     // Sec 4.1.2 / 4.1.5 : 2D Gaussian Ec, diagonal B (streaming or anisotropic diffusion)
#define CR_TEST_CIRCLE_2D         3     // Sec 4.1.5 : 2D ring Ec, circular B         (anisotropic diffusion)
#define CR_TEST_BOTTLENECK_1D     4     // Sec 4.1.3 : 1D density cloud, CR injected from the -x boundary
#define CR_TEST_WAVE_1D           5     // Sec 4.2.1 : 1D CR-driven wave              (full MHD + CR coupling)
#define CR_TEST_BLAST_2D          6     // Sec 4.2.3 : 2D CR-driven blast wave        (full MHD + CR coupling)


// problem-specific global variables
// =======================================================================================
static int    CR_Streaming_Test;        // test selector (one of the CR_TEST_* macros above)
static int    CR_Streaming_Dir;          // streaming/diffusion axis for the 1D tests (0/1/2 = x/y/z)
static double CR_Streaming_FlowV;        // uniform background flow velocity along CR_Streaming_Dir
                                         // (e.g. the moving-fluid diffusion test, Sec 4.1.4)
// =======================================================================================


// function prototypes shared within this file
#if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef COSMIC_RAY
   Aux_Error( ERROR_INFO, "COSMIC_RAY must be enabled !!\n" );
#  endif // #ifndef COSMIC_RAY

#  if ( EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY when enable COSMIC_RAY !!\n" );
#  endif

#  ifdef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be disabled !!\n" );
#  endif

#  ifndef CR_STREAMING
   Aux_Error( ERROR_INFO, "CR_STREAMING must be enabled !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled (CR streaming/diffusion is field-aligned) !!\n" );
#  endif


// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//                3. Must call EoS_Init() before calling any other EoS routine
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CR_Streaming_Test",  &CR_Streaming_Test,      0,           0,                6                 );
   ReadPara->Add( "CR_Streaming_Dir",   &CR_Streaming_Dir,       0,           0,                2                 );
   ReadPara->Add( "CR_Streaming_FlowV", &CR_Streaming_FlowV,     0.0,         NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   const bool is_1D = ( CR_Streaming_Test == CR_TEST_TRIANGULAR_1D  ||
                        CR_Streaming_Test == CR_TEST_GAUSSIAN_1D    ||
                        CR_Streaming_Test == CR_TEST_BOTTLENECK_1D  ||
                        CR_Streaming_Test == CR_TEST_WAVE_1D          );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// for the 1D tests, output a line along the streaming direction
   if ( is_1D )
   {
      const int target = ( CR_Streaming_Dir == 0 ) ? OUTPUT_X : ( CR_Streaming_Dir == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      if ( OPT__OUTPUT_PART != target )
      {
         OPT__OUTPUT_PART = target;
         PRINT_RESET_PARA( OPT__OUTPUT_PART, FORMAT_INT, "" );
      }
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  CR_Streaming_Test     = %d\n",     CR_Streaming_Test  );
      Aux_Message( stdout, "  CR_Streaming_Dir      = %d\n",     CR_Streaming_Dir   );
      Aux_Message( stdout, "  CR_Streaming_FlowV    = %14.7e\n", CR_Streaming_FlowV );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CR_Bottleneck_Dens
// Description :  Density profile of the cold cloud used in the bottleneck test (Sec 4.1.3, Eq. 24)
//
// Note        :  1. Uses ABSOLUTE coordinates (the cloud is centered at x0 = 200 in the domain x in [0,1000])
//                2. rho = rho_h + (rho_c - rho_h)*[1+tanh((x-x0)/dx0)]*[1+tanh((x0-x)/dx0)]
//                   --> a dense cold cloud (rho_c) embedded in a diffuse hot background (rho_h)
//
// Parameter   :  x : Position along the streaming direction (absolute coordinate)
//
// Return      :  gas density
//-------------------------------------------------------------------------------------------------------
static double CR_Bottleneck_Dens( const double x )
{
   const double rho_c = 1.0, rho_h = 0.1, x0 = 200.0, dx0 = 25.0;
   const double t1 = 1.0 + std::tanh( (x  - x0)/dx0 );
   const double t2 = 1.0 + std::tanh( (x0 - x )/dx0 );
   return rho_h + (rho_c - rho_h)*t1*t2;
} // FUNCTION : CR_Bottleneck_Dens



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
//                5. The ADV_* fields (ADV_SIGMA/ADV_VX/ADV_VY/ADV_VZ) are recomputed by CR_UpdateOpacity()
//                   at the start of every fluid solve, so the IC values set here are overwritten before use
//                   --> we therefore set safe no-streaming defaults
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double xc = amr->BoxCenter[0];
   const double yc = amr->BoxCenter[1];

   double Dens = 1.0, vx = 0.0, vy = 0.0, vz = 0.0, Pgas = 1.0, cr_E = 0.0;

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      {
//       Sec 4.1.1: Ec = 2 - |r| for |r| < 1, else 1   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         const double d = std::fabs( r );
         cr_E = ( d < 1.0 ) ? (2.0 - d) : 1.0;
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_1D :
      {
//       Sec 4.1.1 / 4.1.4: Ec = exp(-40 r^2)   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         cr_E = std::exp( -40.0 * r*r );
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_2D :
      {
//       Sec 4.1.2 / 4.1.5: Ec = exp(-40 (dx^2 + dy^2)),  B along the diagonal (set in SetBFieldIC)
         const double dx = x - xc, dy = y - yc;
         cr_E = std::exp( -40.0 * (dx*dx + dy*dy) );
         break;
      }

      case CR_TEST_CIRCLE_2D :
      {
//       Sec 4.1.5: a ring of Ec; B is circular (set in SetBFieldIC)
//       Ec = 12 in 0.5 < r < 0.7 and |phi| < pi/12 (consistent with the analytic solution Eq. 28), else 10
         const double dx = x - xc, dy = y - yc;
         const double r   = std::sqrt( dx*dx + dy*dy );
         const double phi = std::atan2( dy, dx );
         cr_E = ( r > 0.5  &&  r < 0.7  &&  std::fabs(phi) < M_PI/12.0 ) ? 12.0 : 10.0;
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      {
//       Sec 4.1.3: a cold dense cloud; CRs are injected from the -x boundary (see BottleneckBC)
//       --> the simulation domain starts essentially empty of CRs
         Dens = CR_Bottleneck_Dens( x );        // absolute coordinate (cloud at x0 = 200)
         cr_E = 1.0e-6;
         break;
      }

      case CR_TEST_WAVE_1D :
      {
//       Sec 4.2.1: Ec = 20 + 10 sin(pi*(x-xc)); uniform gas; full MHD + CR coupling (CR_SOURCE=1)
         cr_E = 20.0 + 10.0*std::sin( M_PI*(x - xc) );
         break;
      }

      case CR_TEST_BLAST_2D :
      {
//       Sec 4.2.3: Ec = 100 inside r < 0.02, else 0.1; uniform background (rho=1, Eint=2.5)
         const double dx = x - xc, dy = y - yc;
         const double r  = std::sqrt( dx*dx + dy*dy );
         cr_E = ( r < 0.02 ) ? 100.0 : 0.1;
         Pgas = (GAMMA - 1.0)*2.5;               // background internal energy = 2.5
         break;
      }

      default :
         Aux_Error( ERROR_INFO, "unsupported CR_Streaming_Test (%d) !!\n", CR_Streaming_Test );
   } // switch ( CR_Streaming_Test )


// momentum
   const double MomX = Dens*vx;
   const double MomY = Dens*vy;
   const double MomZ = Dens*vz;

// gas + (placeholder) CR pressure for the COSMIC_RAY EoS
// --> CRAY is the EoS cosmic-ray field and is decoupled from the two-moment CR_E used by the
//     streaming module; the gas-CR back-reaction is handled by the streaming source term via CR_E/CR_F
   const double P_cr = 1.0;
   const double Pres = Pgas + P_cr;
   fluid[CRAY] = P_cr / (GAMMA_CR - 1.0);

   const double Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );
   const double Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );   // do NOT include magnetic energy here

// two-moment CR fields (Ec and the three flux components; initial flux = 0 as in the paper)
   fluid[CR_E ] = cr_E;
   fluid[CR_F1] = 0.0;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

// streaming opacity/velocity: recomputed by CR_UpdateOpacity() each step --> safe no-streaming defaults
   fluid[ADV_SIGMA] = CR_MAX_OPACITY;
   fluid[ADV_VX   ] = 0.0;
   fluid[ADV_VY   ] = 0.0;
   fluid[ADV_VZ   ] = 0.0;

// gas fields
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                2. The Alfven speed used by the streaming module is v_A = |B|/sqrt(rho); the field
//                   geometry also sets the direction of field-aligned streaming/diffusion
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      case CR_TEST_GAUSSIAN_1D :
//       uniform field along the streaming direction, |B| = 1  --> v_A = 1
         magnetic[MAGX] = ( CR_Streaming_Dir==0 ) ? 1.0 : 0.0;
         magnetic[MAGY] = ( CR_Streaming_Dir==1 ) ? 1.0 : 0.0;
         magnetic[MAGZ] = ( CR_Streaming_Dir==2 ) ? 1.0 : 0.0;
         break;

      case CR_TEST_GAUSSIAN_2D :
//       uniform field along the x-y diagonal, |B| = 1  --> v_A = 1   (Sec 4.1.2 / 4.1.5)
         magnetic[MAGX] = 1.0/std::sqrt(2.0);
         magnetic[MAGY] = 1.0/std::sqrt(2.0);
         magnetic[MAGZ] = 0.0;
         break;

      case CR_TEST_CIRCLE_2D :
      {
//       circular field B = ( -(y-yc), (x-xc) ) / r, |B| = 1   (Sec 4.1.5)
//       --> analytically divergence-free; corresponds to the vector potential Az = -r
         const double dx = x - amr->BoxCenter[0];
         const double dy = y - amr->BoxCenter[1];
         const double r  = std::sqrt( dx*dx + dy*dy );
         if ( r > TINY_NUMBER ) {
            magnetic[MAGX] = -dy/r;
            magnetic[MAGY] =  dx/r;
         } else {
            magnetic[MAGX] = 0.0;
            magnetic[MAGY] = 0.0;
         }
         magnetic[MAGZ] = 0.0;
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      case CR_TEST_WAVE_1D :
      case CR_TEST_BLAST_2D :
      default :
//       uniform field along x, |B| = 1
         magnetic[MAGX] = 1.0;
         magnetic[MAGY] = 0.0;
         magnetic[MAGZ] = 0.0;
         break;
   } // switch ( CR_Streaming_Test )

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  BottleneckBC
// Description :  User boundary condition for the bottleneck test (Sec 4.1.3): inject CRs from the
//                -x boundary by fixing Ec = 3 while keeping the background gas profile
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//                2. Only the -x face is set to the user BC (OPT__BC_FLU_XM = 4); the +x face uses outflow
//                3. The paper sets the boundary CR flux to the (sign-flipped) value of the last active
//                   zone ("reflecting"); here we use Fc = 0 at the ghost, which is sufficient to inject
//                   CRs from the boundary and drive the steady-state bottleneck profile
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BottleneckBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                   const int GhostSize, const int idx[], const double pos[], const double Time,
                   const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// start from the background IC (gas density profile, Ec = 1e-6, Fc = 0, ADV defaults)
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

// then fix the injected CR energy density and zero the CR flux at the boundary
   fluid[CR_E ] = 3.0;
   fluid[CR_F1] = 0.0;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

} // FUNCTION : BottleneckBC
#endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_Streaming
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_Streaming()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
#  endif

// the bottleneck test injects CRs through a user boundary condition on the -x face
   if ( CR_Streaming_Test == CR_TEST_BOTTLENECK_1D )
   {
      BC_User_Ptr                    = BottleneckBC;
//    with MHD on, a user fluid BC also requires a matching magnetic-field user BC; the field is
//    uniform (Bx = 1) in this test, so SetBFieldIC already provides the correct ghost-zone values
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }
#  endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_Streaming
