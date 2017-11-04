/*!
 * \file braid_wrapper.cpp
 * \brief Functions for XBraid integration.
 *
 * \author S. Guenther
 *
 */

//#include <braid.hpp>
#include <../include/util.hpp>
#include <../include/braid_structure.hpp>
#include "../../Common/include/ad_structure.hpp"


//BraidTape_t* braidTape;

//void setupTapeData(){
//   /* Create the braid tape */
//   braidTape = new BraidTape_t();
//}


/* Make a copy of a Vector. Used for primal taping. */
braid_Vector deep_copy( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(nPoint, nVar);

    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Copy the values from u to v */
        /* TODO: Set the values with the SetSolution(...), Set_solution_time_n(...) and Set_Solution_timen1(...) */
        for (int iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
        }
    }
    return v;
}


int my_Step( braid_App        app,
             braid_Vector     ustop,
             braid_Vector     fstop,
             braid_Vector     u,
             braid_StepStatus status ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    int nDim   = app->geometry_container[ZONE_0][MESH_0]->GetnDim();

    /* Grab status of current time step from xBraid */
    double tstart;
    double tstop;
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
////    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

    /* Declare and allocate intermediate casting variables */
    su2double* cast_n  = new su2double[nVar];
    su2double* cast_n1 = new su2double[nVar];


    /* Trick SU2 with xBraid's DeltaT / 2 */
    double deltat = ( tstop - tstart ) / 2.0;
    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );

    /* Trick the su2 solver with the correct state vector (Solution, Solution_time_n and Solution_time_n1*/
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      for (int iVar = 0; iVar < nVar; iVar++){
        cast_n[iVar]  = u->Solution->time_n[iPoint][iVar];
        cast_n1[iVar] = u->Solution->time_n1[iPoint][iVar];
      }
      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast_n);
      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(cast_n);
      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(cast_n1);
    }

    /* Trick SU2 with the correct iExtIter = (t - t0)/dt  -1 */
    int iExtIter = (int) round( ( tstart + deltat - app->initialstart ) / app->initialDT) -1 ;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Print information output */
     if (app->su2rank == MASTER_NODE)
         cout << "rank " << app->braidrank << ": " << deltat << "-step from " << tstart << " to " << tstop << endl;

    /* Take the first time step to tstart + deltat */
    app->driver->Run();

    /* Update the Solution_n and solution_n1 for dual time stepping strategy */
    app->driver->Update();


    /* Grab drag and lift coefficient from SU2's master node. */
    if (app->su2rank == MASTER_NODE){
        u->Solution->Total_CD_n1 = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CD());
        u->Solution->Total_CL_n1 = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CL());
      }


     /* Trick SU2 with the next iExtIter */
     iExtIter++;
     app->config_container[ZONE_0]->SetExtIter(iExtIter);

     /* Take the next time step to tstart + 2*deltat = tstop */
     app->driver->Run();

     /* Update the Solution_n and solution_n1 for dual time stepping strategy */
     app->driver->Update();

     /* Grab the history values from SU2's master node. */
     if (app->su2rank == MASTER_NODE){

       u->Solution->Total_CD_n = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CD());
       u->Solution->Total_CL_n = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CL());

     }


     /* Grab the solution vectors from su2 for both time steps */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
            u->Solution->time_n1[iPoint][iVar] = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
        }
    }

    /* Free memory of the intermediate casting variables */
    delete [] cast_n;
    delete [] cast_n1;

    /* Tell XBraid no refinement */
    braid_StepStatusSetRFactor(status, 1);

    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    double Density_Inf   = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetDensity_FreeStreamND());
    su2double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
    double Energy_Inf    = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetEnergy_FreeStreamND());
    double Pressure_Inf  = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetPressure_FreeStreamND());
    bool compressible = (app->config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (app->config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);

    /* Allocate memory for the primal braid vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(nPoint, nVar);

    /*--- Initialize the solution to the far-field state everywhere. ---*/

    for (int iPoint = 0; iPoint < nPoint; iPoint++){

        /* Initialize the solution with the freestream values */
        if (compressible) {
			u->Solution->time_n[iPoint][0]  = Density_Inf;
			u->Solution->time_n1[iPoint][0] = Density_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution->time_n[iPoint][iDim+1]  = Density_Inf*SU2_TYPE::GetValue(Velocity_Inf[iDim]);
				u->Solution->time_n1[iPoint][iDim+1] = Density_Inf*SU2_TYPE::GetValue(Velocity_Inf[iDim]);
			}
			u->Solution->time_n[iPoint][nVar-1]  = Density_Inf*Energy_Inf;
            u->Solution->time_n1[iPoint][nVar-1] = Density_Inf*Energy_Inf;
		}
        if (incompressible) {
			u->Solution->time_n[iPoint][0]  = Pressure_Inf;
			u->Solution->time_n1[iPoint][0] = Pressure_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution->time_n[iPoint][iDim+1]  = SU2_TYPE::GetValue(Velocity_Inf[iDim])*Density_Inf;
				u->Solution->time_n1[iPoint][iDim+1] = SU2_TYPE::GetValue(Velocity_Inf[iDim])*Density_Inf;
			}
		}
	}

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}

int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr ){

    /* Grab variables from the app */
    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nDim        = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    CConfig *config = app->config_container[ZONE_0];

    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(nPoint, nVar);

    /* Copy the values from u to v */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
        }
    }

    /* Set the pointer */
    *v_ptr = v;

    return 0;
}

int my_Free( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    /* Delete the primal solution lists (calls destructor of TwoStepSolution) */
    delete u->Solution;

    /* Delete the braid_Vector */
    delete u;

    return 0;
}

int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute the sum y = alpha x + beta y at time n and time n-1 */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
                                              + beta  * y->Solution->time_n[iPoint][iVar];
            y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
                                              + beta  * y->Solution->time_n1[iPoint][iVar];
        }
    }

    return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute l2norm of the solution list at time n and n1 */
    double norm = 0.0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->Solution->time_n[iPoint][iVar], 2);
            norm += pow(u->Solution->time_n1[iPoint][iVar], 2);
        }
    }

    /* Set the pointer */
    *norm_ptr = sqrt(norm);

    return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){


    /* --- Add drag and lift to objective function. ---*/
    app->Total_CD_avg += u->Solution->Total_CD_n + u->Solution->Total_CD_n1;
    app->Total_CL_avg += u->Solution->Total_CL_n + u->Solution->Total_CL_n1;



    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /*--- Write solution output files only if XBraid has finished ---*/
    if (app->done) {

      /* Retrieve xBraid information from status object */
      double t;
      braid_AccessStatusGetT(astatus, &t);
  
      /* Allocate memory for the casting variable */
      su2double* cast = new su2double[nVar];
  
      /* Compute the time step iExtIter = (t - t0)/dt -1 which is used for naming the restart file */
      int iExtIter = (int) round( ( t - app->initialstart ) / app->initialDT) -1 ;
      app->config_container[ZONE_0]->SetExtIter(iExtIter);
  
  
      /* Only continue if iExtIter > 0 !! Otherwise xbraid tries to write at timestep -1*/
      if (iExtIter>0){
  
        /* --- Write Solution_time_n to restart file ---*/
  
        /* Trick SU2 with the current solution for output (SU2 writes CVariable::Solution, not _time_n!) */
        for (int iPoint = 0; iPoint < nPoint; iPoint++){
          for (int iVar = 0.0; iVar < nVar; iVar++){
            cast[iVar] = u->Solution->time_n[iPoint][iVar];
          }
          app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
        }
  
        /* Compute the primitive Variables from the conservative ones */
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);
  
  //    /* Call the SU2 output routine */
        if (app->su2rank==MASTER_NODE)
            cout << "rank_t " << app->braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
        app->output->SetResult_Files_Parallel(app->solver_container, app->geometry_container,
                                   app->config_container, iExtIter, 1);
  
        /* Write history values at time n to the app stream */
        /* NOT WORKING RIGHT NOW....... */
        //if (app->su2rank == MASTER_NODE){
        //    *app->history_stream << iExtIter << " " << u->Solution->Total_CLift_n
        //                         << " " << u->Solution->Total_CDrag_n
        //                         << " " << u->Solution->Total_CSideForce_n
        //                         << " " << u->Solution->Total_CMx_n
        //                         << " " << u->Solution->Total_CMy_n
        //                         << " " << u->Solution->Total_CMz_n
        //                         << " " << u->Solution->Total_CFx_n
        //                         << " " << u->Solution->Total_CFy_n
        //                         << " " << u->Solution->Total_CFz_n
        //                         << " " << u->Solution->Total_CEff_n;
        //    //      for (int iVar = 0; iVar < nVar; iVar++){
        //    //        *app->history_stream << " " << u->residual_flow_n[iVar];
        //    //      }
        //    *app->history_stream << " " << u->residual_dens_n;
        //    *app->history_stream << "\n";
        //}
  
  
        /* --- Write Solution_time_n1 to restart file ---*/
  
  
        /* Trick SU2 with the correct iExtIter = iExtIter - 1 */
        iExtIter--;
        app->config_container[ZONE_0]->SetExtIter(iExtIter);
  
        /* Trick SU2 with the current solution for output */
        for (int iPoint = 0; iPoint < nPoint; iPoint++){
          for (int iVar = 0.0; iVar < nVar; iVar++){
            cast[iVar] = u->Solution->time_n1[iPoint][iVar];
          }
          app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
        }
  
        /* Compute the primitive Variables from the conservative ones */
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);
  
        /* Call the SU2 output routine */
        if (app->su2rank==MASTER_NODE) cout << "rank_t " << app->braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
        app->output->SetResult_Files_Parallel(app->solver_container, app->geometry_container,
                                         app->config_container, iExtIter, 1);
  
        /* Write history values at time n-1 to the app stream */
        /* NOT WORKING RIGHT NOW....... */
        //if (app->su2rank == MASTER_NODE){
        //    *app->history_stream << iExtIter << " " << u->Solution->Total_CLift_n1
        //                         << " " << u->Solution->Total_CDrag_n1
        //                         << " " << u->Solution->Total_CSideForce_n1
        //                         << " " << u->Solution->Total_CMx_n1
        //                         << " " << u->Solution->Total_CMy_n1
        //                         << " " << u->Solution->Total_CMz_n1
        //                         << " " << u->Solution->Total_CFx_n1
        //                         << " " << u->Solution->Total_CFy_n1
        //                         << " " << u->Solution->Total_CFz_n1
        //                         << " " << u->Solution->Total_CEff_n1;
        //    //      for (int iVar = 0; iVar < nVar; iVar++){
        //    //        *app->history_stream << " " << u->residual_flow_n1[iVar];
        //    //      }
        //    *app->history_stream << " " << u->residual_dens_n1;
        //    *app->history_stream << "\n";
        //}
  
  
  
        /* Free memory for the intermediat casting variable */
        delete [] cast;
  
      }
    }

    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute size of buffer */
    *size_ptr = 2.0 * nPoint * nVar * sizeof(double);

    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Pack the buffer with current and previous time */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Write Solution at current time to the buffer */
          dbuffer[ibuffer] = u->Solution->time_n[iPoint][iVar];
          ibuffer++;
          /* Write Solution at previous time to the buffer */
          dbuffer[ibuffer] = u->Solution->time_n1[iPoint][iVar];
          ibuffer++;
        }
    }

    /* Set the size of the buffer */
    int size = 2.0 * nPoint * nVar * sizeof(double);
    braid_BufferStatusSetSize( bstatus, size );


    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Get the buffer */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;

    /* Allocate memory for the new braid Vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(nPoint, nVar);

    /* Unpack the buffer and write solution at current and previous time */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          u->Solution->time_n[iPoint][iVar] = dbuffer[ibuffer];
          ibuffer++;
          u->Solution->time_n1[iPoint][iVar] = dbuffer[ibuffer];
          ibuffer++;
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}
