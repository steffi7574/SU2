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

std::vector<void*>* tape;

//void setupTapeData() {
//    tape = new std::vector<void*>();
//    tape->resize(100);
//}

int my_Step( braid_App        app,
             braid_Vector     ustop,
             braid_Vector     fstop,
             braid_Vector     u,
             braid_StepStatus status ){
//int my_Phi( braid_App app, braid_Vector u, braid_PhiStatus status ){

    /* Grab rank of the current SU2 processor */


    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Grab status of current time step from xBraid */
    su2double tstart;
    su2double tstop;
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

    /* Trick SU2 with xBraid's DeltaT / 2 */
    su2double deltat = ( tstop - tstart ) / 2.0;
    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );

    /* Trick the su2 solver with the correct state vector (Solution, Solution_time_n and Solution_time_n1*/
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution_time_n[iPoint]);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(u->Solution_time_n[iPoint]);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(u->Solution_time_n1[iPoint]);
    }

    /* Trick SU2 with the correct iExtIter = (t - t0)/dt  -1 */
    int iExtIter = (int) round( ( tstart + deltat - app->initialstart ) / app->initialDT) ;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Print information output */
    if (app->su2rank == MASTER_NODE) cout<<format(" %d: two %1.3f-steps from %1.4f to %1.4f\n", app->braidrank, deltat, tstart, tstop);

    /* Take the first time step to tstart + deltat */
//    app->driver->Run(app->iteration_container, app->output, app->integration_container,
//                   app->geometry_container, app->solver_container, app->numerics_container,
//                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//                   app->interpolator_container, app->transfer_container);

    /* Grab the history values from SU2's master node. */
    if (app->su2rank == MASTER_NODE){
//      /* Grab the flow coefficient values for that time step and store it at time n-1 */
      u->Total_CLift_n1       = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift();
      u->Total_CDrag_n1       = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag();
      u->Total_CSideForce_n1  = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CSideForce();
      u->Total_CEff_n1        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CEff();
      u->Total_CMx_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMx();
      u->Total_CMy_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMy();
      u->Total_CMz_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMz();
      u->Total_CFx_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFx();
      u->Total_CFy_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFy();
      u->Total_CFz_n1         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFz();
//      /* Grab the flow residual at that time stap and store it at time n-1 */
//      for (int iVar = 0; iVar < nVar; iVar++){
//          cout << "HI " << iVar <<"\n";
//         u->residual_flow_n1[iVar] = 3.2;
//         cout << "HI" << u->residual_flow_n1[iVar] <<"\n";
//         u->residual_flow_n1[iVar] = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(iVar);
//      }
      u->residual_dens_n1 = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0);
    }


    /* Trick SU2 with the next iExtIter */
    iExtIter++;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Take the next time step to tstart + 2*deltat = tstop */
//    app->driver->Run(app->iteration_container, app->output, app->integration_container,
//                   app->geometry_container, app->solver_container, app->numerics_container,
//                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//                   app->interpolator_container, app->transfer_container);

    /* Grab the history values from SU2's master node. */
    if (app->su2rank == MASTER_NODE){
      /* Grab the flow coefficient values for that time step and store it at time n-1 */
      u->Total_CLift_n       = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift();
      u->Total_CDrag_n       = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag();
      u->Total_CSideForce_n  = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CSideForce();
      u->Total_CEff_n        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CEff();
      u->Total_CMx_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMx();
      u->Total_CMy_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMy();
      u->Total_CMz_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMz();
      u->Total_CFx_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFx();
      u->Total_CFy_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFy();
      u->Total_CFz_n         = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFz();
      /* Grab the flow residual at that time stap and store it at time n-1 */
//      for (int iVar = 0; iVar < nVar; iVar++){
//         u->residual_flow_n[iVar] = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(iVar);
//      }
      u->residual_dens_n = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0);
    }


   /* Grab the solution vectors from su2 for both time steps */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            u->Solution_time_n[iPoint][iVar]  = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar];
            u->Solution_time_n1[iPoint][iVar] = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar];
        }
    }

    /* Tell XBraid no refinement */
    braid_StepStatusSetRFactor(status, 1);


    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    su2double Density_Inf   = app->config_container[ZONE_0]->GetDensity_FreeStreamND();
    su2double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
    su2double Energy_Inf    = app->config_container[ZONE_0]->GetEnergy_FreeStreamND();
    su2double Pressure_Inf  = app->config_container[ZONE_0]->GetPressure_FreeStreamND();
    bool compressible = (app->config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (app->config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
    bool freesurface = (app->config_container[ZONE_0]->GetKind_Regime() == FREESURFACE);

    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      /* Print information output */
      if (app->su2rank == MASTER_NODE) cout << format("My_Init at time %1.4f\n", t);
    }

    /* Allocate memory */
    my_Vector* u;
    u = new my_Vector;
    u->Solution_time_n  = new su2double*[nPoint];
    u->Solution_time_n1 = new su2double*[nPoint];

    for (int iPoint = 0; iPoint < nPoint; iPoint++){

        /* Allocate memory */
        u->Solution_time_n[iPoint]  = new su2double[nVar];
        u->Solution_time_n1[iPoint] = new su2double[nVar];

        /* Initialize the solution with the freestream values */
        if (compressible) {
			u->Solution_time_n[iPoint][0]  = Density_Inf;
			u->Solution_time_n1[iPoint][0] = Density_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution_time_n[iPoint][iDim+1]  = Density_Inf*Velocity_Inf[iDim];
				u->Solution_time_n1[iPoint][iDim+1] = Density_Inf*Velocity_Inf[iDim];
			}
			u->Solution_time_n[iPoint][nVar-1]  = Density_Inf*Energy_Inf;
            u->Solution_time_n1[iPoint][nVar-1] = Density_Inf*Energy_Inf;
		}
        if (incompressible || freesurface) {
			u->Solution_time_n[iPoint][0]  = Pressure_Inf;
			u->Solution_time_n1[iPoint][0] = Pressure_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution_time_n[iPoint][iDim+1]  = Velocity_Inf[iDim]*Density_Inf;
				u->Solution_time_n1[iPoint][iDim+1] = Velocity_Inf[iDim]*Density_Inf;
			}
		}
	}


//    /* Initialize Flow residual */
//    u->residual_flow_n  = new su2double[nVar];
//    u->residual_flow_n1 = new su2double[nVar];
//    for (int iVar = 0; iVar < nVar; iVar++){
//       u->residual_flow_n[iVar]  = 0.0;
//       u->residual_flow_n1[iVar] = 0.0;
//    }

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

    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      /* Print information output */
      if (app->su2rank == MASTER_NODE) cout << format("My_Clone\n");
    }

    /* Allocate memory for the new copy v */
    my_Vector* v;
    v = new my_Vector;
    v->Solution_time_n = new su2double*[nPoint];
    v->Solution_time_n1 = new su2double*[nPoint];

    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        v->Solution_time_n[iPoint]  = new su2double[nVar];
        v->Solution_time_n1[iPoint] = new su2double[nVar];
        /* Copy the values from u to v */
        for (int iVar = 0; iVar < nVar; iVar++){
            v->Solution_time_n[iPoint][iVar]  = u->Solution_time_n[iPoint][iVar];
            v->Solution_time_n1[iPoint][iVar] = u->Solution_time_n1[iPoint][iVar];
        }
    }

    /* Set the pointer */
    *v_ptr = v;

    return 0;
}

int my_Free( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      /* Print information output */
      if (app->su2rank == MASTER_NODE) cout << format("My_Free\n");
    }

    /* Delete the Solution each point in space. */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        delete [] u->Solution_time_n[iPoint];
        delete [] u->Solution_time_n1[iPoint];
    }

    /* Delete the list of Solutions */
    delete [] u->Solution_time_n;
    delete [] u->Solution_time_n1;

//    /* Delete the flow residual */
//    delete [] u->residual_flow_n;
//    delete [] u->residual_flow_n1;

    /* Delete braid vector */
    delete u;

    return 0;
}

int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      /* Print information output */
      if (app->su2rank == MASTER_NODE) cout << format("My_Sum\n");
    }

    /* Loop over all points */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Loop over all variables */
        for (int iVar = 0; iVar < nVar; iVar++){
            /* Compute the sum y = alpha x + beta y at time n and time n-1 */
            y->Solution_time_n[iPoint][iVar]  = alpha * x->Solution_time_n[iPoint][iVar]
                                              + beta  * y->Solution_time_n[iPoint][iVar];
            y->Solution_time_n1[iPoint][iVar] = alpha * x->Solution_time_n1[iPoint][iVar]
                                              + beta  * y->Solution_time_n1[iPoint][iVar];
        }
    }

    return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      /* Print information output */
      if (app->su2rank == MASTER_NODE) cout << format("My_SpatialNorm\n");
    }

    /* Compute l2norm of the solution list at time n and n1 */
    su2double norm = 0.0;
    /* Loop over all points */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Loop over all variables */
        for (int iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->Solution_time_n[iPoint][iVar], 2);
            norm += pow(u->Solution_time_n1[iPoint][iVar], 2);
        }
    }

    /* Set the pointer */
    *norm_ptr = sqrt(norm);

    return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Retrieve xBraid time information from status object */
    double t;
    braid_AccessStatusGetT(astatus, &t);

    /* Print information output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_Access at t = %1.4f\n", t);
    }

    /* --- Write Solution_time_n to restart file ---*/


    /* Trick SU2 with the correct iExtIter = (t - t0)/dt which is used for naming the restart file */
    int iExtIter = (int) round( ( t - app->initialstart ) / app->initialDT) ;
//    iExtIter = iExtIter + 10000*braidrank;   /* Make the identifier unique for each braid processor */
    app->config_container[ZONE_0]->SetExtIter(iExtIter);
    /* Check if xBraid tries to write to negative iExtIter */
    if (iExtIter < 0) {
      if (app->su2rank==MASTER_NODE) cout << "rank_t " << app->braidrank << " tries to write to iExtIter -1 -> Early Exit.\n" << endl;
      return 0;
    }

    /* Trick SU2 with the current solution for output (SU2 writes CVariable::Solution, not _time_n!) */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution_time_n[iPoint]);
    }

    /* Compute the primitive Variables from the conservative ones */
    app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

//    /* Call the SU2 output routine */
//    if (app->su2rank==MASTER_NODE) cout << "rank_t " << braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
//    app->output->SetResult_Files(app->solver_container, app->geometry_container,
//                                 app->config_container, iExtIter, 1);

    /* Write history values at time n to the app stream */
    if (app->su2rank == MASTER_NODE){
      *app->history_stream << iExtIter << " " << u->Total_CLift_n
                                       << " " << u->Total_CDrag_n
                                       << " " << u->Total_CSideForce_n
                                       << " " << u->Total_CMx_n
                                       << " " << u->Total_CMy_n
                                       << " " << u->Total_CMz_n
                                       << " " << u->Total_CFx_n
                                       << " " << u->Total_CFy_n
                                       << " " << u->Total_CFz_n
                                       << " " << u->Total_CEff_n;
//      for (int iVar = 0; iVar < nVar; iVar++){
//        *app->history_stream << " " << u->residual_flow_n[iVar];
//      }
      *app->history_stream << " " << u->residual_dens_n;
      *app->history_stream << "\n";
    }

    /* --- Write Solution_time_n1 to restart file ---*/


    /* Trick SU2 with the correct iExtIter = iExtIter - 1 */
    iExtIter--;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);
    /* Check if xBraid tries to write to negative iExtIter */
    if (iExtIter < 0) {
      if (app->su2rank==MASTER_NODE) cout << "rank_t " << app->braidrank << " tries to write to iExtIter -1 -> Early Exit.\n" << endl;
      return 0;
    }

    /* Trick SU2 with the current solution for output */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution_time_n1[iPoint]);
    }

    /* Compute the primitive Variables from the conservative ones */
    app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

//    /* Call the SU2 output routine */
//    app->output->SetResult_Files(app->solver_container, app->geometry_container,
//                                 app->config_container, iExtIter, 1);


    /* Write history values at time n-1 to the app stream */
    if (app->su2rank == MASTER_NODE){
      *app->history_stream << iExtIter << " " << u->Total_CLift_n1
                                       << " " << u->Total_CDrag_n1
                                       << " " << u->Total_CSideForce_n1
                                       << " " << u->Total_CMx_n1
                                       << " " << u->Total_CMy_n1
                                       << " " << u->Total_CMz_n1
                                       << " " << u->Total_CFx_n1
                                       << " " << u->Total_CFy_n1
                                       << " " << u->Total_CFz_n1
                                       << " " << u->Total_CEff_n1;
//      for (int iVar = 0; iVar < nVar; iVar++){
//        *app->history_stream << " " << u->residual_flow_n1[iVar];
//      }
      *app->history_stream << " " << u->residual_dens_n1;
      *app->history_stream << "\n";
    }

    // sum up the drag value into app->sum



    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute size of buffer */
    *size_ptr = 2.0 * nPoint * nVar * sizeof(su2double);

    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Print information output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_BufPack\n");
    }

    /* Pack the buffer with current and previous time */
    su2double *dbuffer = (su2double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Write Solution at current time to the buffer */
          dbuffer[ibuffer] = u->Solution_time_n[iPoint][iVar];
          ibuffer++;
          /* Write Solution at previous time to the buffer */
          dbuffer[ibuffer] = u->Solution_time_n1[iPoint][iVar];
          ibuffer++;
        }
    }

    /* Set the size of the buffer */
    int size = 2.0 * nPoint * nVar * sizeof(su2double);
    braid_BufferStatusSetSize( bstatus, size );


    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Print information output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_BufUnpack\n");
    }

    /* Get the buffer */
    su2double *dbuffer = (su2double*)buffer;
    int ibuffer = 0;

    /* Allocate memory for the new braid Vector */
    my_Vector* u;
    u = new my_Vector;
    u->Solution_time_n  = new su2double*[nPoint];
    u->Solution_time_n1 = new su2double*[nPoint];

    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        u->Solution_time_n[iPoint]  = new su2double[nVar];
        u->Solution_time_n1[iPoint] = new su2double[nVar];

        /* Unpack the buffer and write solution to current and previous time */
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Unpack Solution at current time from the buffer */
          u->Solution_time_n[iPoint][iVar] = dbuffer[ibuffer];
          ibuffer++;
          /* Unpack Solution at previous time from the buffer */
          u->Solution_time_n1[iPoint][iVar] = dbuffer[ibuffer];
          ibuffer++;
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}
