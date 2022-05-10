#include "GAMER.h"

#if ( defined PARTICLE && defined TRACER )

void Par_MapMesh2Particles ( const double EdgeL[3], const double EdgeR[3], const double _dh,
                             const int AttrSize3D, const real *Attr, const int NPar,
                             real *InterpParPos[3], const real ParType[],
                             const long ParList[], bool useTracers, real ParAttr[],
                             const bool CorrectVelocity );


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateTracerParticle
// Description :  Update tracer particle position and velocity at the target level
//
// Note        :  1. Does not take into account the "periodic B.C." when updating particles
//                   --> After update, the particle position may lie outside the simulation box
//                   --> It will be corrected by Par_PassParticle2Sibling()
//                2. Particle time may not be synchronized (so different particles may have different time).
//                   --> For example, particles just cross from coarse (lv) to fine (lv+1) grids may have time greater than
//                       other particles at lv+1. Also, particles just cross from fine (lv) to coarse (lv-1) grids may have
//                       time less than other particles at lv-1.
//                   --> Only update particles with time < TimeNew
//                3. The mapOnly mode maps the velocity to the particles, but does not advance
//                   their positions
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                TimeOld      : Physical time before update
//                mapOnly      : If true, only map the velocity to the particles
//                               and do not advance them
//-------------------------------------------------------------------------------------------------------
void Par_UpdateTracerParticle( const int lv, const double TimeNew, const double TimeOld,
                               bool mapOnly )
{

   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const double dh                = amr->dh[lv];
   const double _dh               = 1.0/dh;

   const int  ParGhost            = amr->Par->GhostSizeTracer;
   const int  VelSize             = PS1 + 2*ParGhost;

   real *ParPos[3]       = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *ParVel[3]       = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   real *ParTime         = amr->Par->Time;
   const real *ParType   = amr->Par->Type;

// check
#  ifdef COMOVING
#  error : ERROR : does not support COMOVING because time-step has not been converted to comoving !!
#  endif

// OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   real *VelX = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group
   real *VelY = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group
   real *VelZ = new real [ 8*CUBE(VelSize) ];    // 8: number of patches per patch group

   bool   GotYou;
   long   ParID;
   real   dt;

// loop over all **real** patch groups
#  pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with target tracer particles
//    --> use patch group as the calculation unit since Prepare_PatchData() only work with patch group
//    --> disadvantage: some patches may not have tracer particles ... (they will be skipped later)
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         if ( amr->patch[0][lv][PID]->NParType[(long)PTYPE_TRACER] > 0 )
         {
            GotYou = true;
         }

         if ( GotYou )  break;
      } // for (int PID=PID0; PID<PID0+8; PID++)


//    nothing to do if there are no target tracer particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the velocity data for the patch group with particles (need NSIDE_26 for ParGhost>0 )

      Prepare_PatchData( lv, TimeNew, VelX, NULL, ParGhost, 1, &PID0, _VELX, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
      Prepare_PatchData( lv, TimeNew, VelY, NULL, ParGhost, 1, &PID0, _VELY, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,DE_Consistency_No );
      Prepare_PatchData( lv, TimeNew, VelZ, NULL, ParGhost, 1, &PID0, _VELZ, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
//       skip patches with no tracer particles
         if ( amr->patch[0][lv][PID]->NParType[(long)PTYPE_TRACER] == 0 ) continue;

         double EdgeL[3], EdgeR[3];

         real** Vel_Temp = new real*[3];
         real** InterpParPos = new real*[3];
         for (int d=0; d<3; d++) {
            Vel_Temp[d] = new real[amr->patch[0][lv][PID]->NPar];
            InterpParPos[d] = new real[amr->patch[0][lv][PID]->NPar];
            EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - dh*ParGhost;
            EdgeR[d] = amr->patch[0][lv][PID]->EdgeR[d] + dh*ParGhost;
         }

         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip massive particles
            if ( ParType[ParID] != PTYPE_TRACER )
               continue;

//          determine time-step
            dt = (real)TimeNew - ParTime[ParID];

            for (int d=0; d<3; d++)
               if ( mapOnly )
                  InterpParPos[d][p] = ParPos[d][ParID];
               else
                  InterpParPos[d][p] = ParPos[d][ParID]+dt*ParVel[d][ParID];

         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)

         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelX+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, true, Vel_Temp[0],
                                amr->Par->TracerVelCorr );
         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelY+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, true, Vel_Temp[1],
                                amr->Par->TracerVelCorr );
         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelZ+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, true, Vel_Temp[2],
                                amr->Par->TracerVelCorr );

//       5. update particles

         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip massive particles
            if ( ParType[ParID] != PTYPE_TRACER )
               continue;

            if ( mapOnly )
            {

               for (int d=0; d<3; d++)
                  ParVel[d][ParID] = Vel_Temp[d][p];

            }

//          5.0 Euler method
            else if ( amr->Par->IntegTracer == TRACER_INTEG_EULER )
            {

               for (int d=0; d<3; d++) {
                  ParPos[d][ParID] = InterpParPos[d][p];
                  ParVel[d][ParID] = Vel_Temp[d][p];
               }
               ParTime[ParID] = TimeNew;

            }

//          5.1 RK2 scheme
            else if ( amr->Par->IntegTracer == TRACER_INTEG_RK2 )
            {

//             determine time-step
               dt      = (real)TimeNew - ParTime[ParID];

               for (int d=0; d<3; d++)
                  InterpParPos[d][p] = ParPos[d][ParID] +
                     0.5*dt*(Vel_Temp[d][p]+ParVel[d][ParID]);

            } // amr->Par->IntegTracer

         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)`

         if ( !mapOnly && amr->Par->IntegTracer == TRACER_INTEG_RK2 )
         {

            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelX+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, true, Vel_Temp[0],
                                   amr->Par->TracerVelCorr );
            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelY+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, true, Vel_Temp[1],
                                   amr->Par->TracerVelCorr );
            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelZ+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, true, Vel_Temp[2],
                                   amr->Par->TracerVelCorr );

            for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++) {

               ParID = amr->patch[0][lv][PID]->ParList[p];

//             skip massive particles
               if ( ParType[ParID] != PTYPE_TRACER )
                  continue;

               for (int d=0; d<3; d++) {
                  ParPos[d][ParID] = InterpParPos[d][p];
                  ParVel[d][ParID] = Vel_Temp[d][p];
               }
               ParTime[ParID] = TimeNew;

            } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)`

         } // if ( amr->Par->IntegTracer == TRACER_INTEG_RK2 )

         for (int d=0; d<3; d++) {
            delete [] Vel_Temp[d];
            delete [] InterpParPos[d];
         }
         delete [] Vel_Temp;
         delete [] InterpParPos;

      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)

// 6. free memory
   delete [] VelX;
   delete [] VelY;
   delete [] VelZ;

   } // end of OpenMP parallel region

} // FUNCTION : Par_UpdateTracerParticle

#endif // #if defined( PARTICLE ) && defined ( TRACER )