#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double CR_v0;                    // velocity on the streaming direction
static int    CR_Streaming_Dir;         // streaming direction
// =======================================================================================




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

/*#  ifdef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be disabled !!\n" );
#  endif*/

#  ifndef CR_STREAMING
   Aux_Error( ERROR_INFO, "CR_STREAMING must be enabled !!\n" );
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
   ReadPara->Add( "CR_v0",             &CR_v0,                  0.0,         0.0,              NoMax_double      );
   ReadPara->Add( "CR_Streaming_Dir",  &CR_Streaming_Dir,       0,           0,                2                 );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


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

   if (  ( CR_Streaming_Dir == 0 && OPT__OUTPUT_PART != OUTPUT_X )  ||
         ( CR_Streaming_Dir == 1 && OPT__OUTPUT_PART != OUTPUT_Y )  ||
         ( CR_Streaming_Dir == 2 && OPT__OUTPUT_PART != OUTPUT_Z )    )
   {
      OPT__OUTPUT_PART = ( CR_Streaming_Dir == 0 ) ? OUTPUT_X : ( CR_Streaming_Dir == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      PRINT_RESET_PARA( OPT__OUTPUT_PART, FORMAT_INT, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  CR_v0                 = %14.7e\n", CR_v0 );
      Aux_Message( stdout, "  CR_Streaming_Dir      = %d\n",     CR_Streaming_Dir      );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot, P_cr, CRay, cr_E, cr_F1, cr_F2, cr_F3;
   double r;

   switch ( CR_Streaming_Dir )
   {
      case 0: r=x; break;
      case 1: r=y; break;
      case 2: r=z; break;
   } // switch ( CR_Streaming_Dir )

   Dens = 1.0;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Pres = 1.0;

   const double d = std::abs( r - amr->BoxCenter[CR_Streaming_Dir]);
   cr_E = (d < 1) ? 2 - d : 1;

   P_cr = 1.0; // placeholder for now

   const double cr_F = CR_v0*4.0*cr_E/(3.0*CR_VMAX) - ((r < amr->BoxCenter[CR_Streaming_Dir]) ? 1.0 : -1.0)/CR_SIGMA;
   cr_F1 = (CR_Streaming_Dir == 0) ? cr_F : 0.0;
   cr_F2 = (CR_Streaming_Dir == 1) ? cr_F : 0.0;
   cr_F3 = (CR_Streaming_Dir == 2) ? cr_F : 0.0;

   const double GAMMA_CR_m1_inv = 1.0 / (GAMMA_CR - 1.0);
   Pres = Pres + P_cr;
   CRay = GAMMA_CR_m1_inv*P_cr;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

// set the output array of passive scaler
   fluid[CRAY] = CRay;

   fluid[CR_E] = cr_E;
   fluid[CR_F1] = cr_F1;
   fluid[CR_F2] = cr_F2;
   fluid[CR_F3] = cr_F3;

// set the output array
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

   const double CR_BField = 1.0;
   magnetic[MAGX] = CR_Streaming_Dir == 0 ? CR_BField : 0.0;
   magnetic[MAGY] = CR_Streaming_Dir == 1 ? CR_BField : 0.0;
   magnetic[MAGZ] = CR_Streaming_Dir == 2 ? CR_BField : 0.0;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
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
#  endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY  &&  defined CR_STREAMING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_Streaming
