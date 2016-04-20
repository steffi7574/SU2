/*!
 * \file braid_wrapper.cpp
 * \brief Functions for XBraid integration.
 *
 * \author S. Guenther
 *
 */

//#include <braid.hpp>
#include <../include/braid_structure.hpp>

int my_Phi( braid_App app, braid_Vector u, braid_PhiStatus status ){

    /* Grab rank of the current SU2 processor */
    int su2rank, su2size;
    int braidrank;
    MPI_Comm_rank(app->comm_x, &su2rank);
    MPI_Comm_size(app->comm_x, &su2size);
    MPI_Comm_rank(app->comm_t, &braidrank);


    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    /* Grab status of current time step from xBraid */
    su2double tstart;
    su2double tstop;
    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

    /* Trick SU2 with xBraid's DeltaT / 2 */
    su2double deltat = ( tstop - tstart ) / 2.0;
    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );

    /* Trick the su2 solver with the right state vector (Solution, Solution_time_n and Solution_time_n1*/
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Get the Solution from braid vector u */
        su2double* uSolution_time_n  = u->node[iPoint]->GetSolution_time_n();
        su2double* uSolution_time_n1 = u->node[iPoint]->GetSolution_time_n1();
        /* Set the solution to the SU2 solver */
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(uSolution_time_n);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(uSolution_time_n);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(uSolution_time_n1);
    }

    /* Trick SU2 with the correct iExtIter = (t - t0)/dt - 1  */
    int iExtIter = (int) round( ( tstart + deltat - app->initialstart ) / app->initialDT) - 1;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Print information output */
    if (su2rank == MASTER_NODE) {
      cout << "rank_t " << braidrank << " performes two " << deltat << "-steps from " << tstart << " to " << tstop << endl;
    }

    /* Take the first time step to tstart + deltat */
    app->driver->Run(app->iteration_container, app->output, app->integration_container,
                   app->geometry_container, app->solver_container, app->numerics_container,
                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
                   app->interpolator_container, app->transfer_container);

    /* Trick SU2 with the next iExtIter */
    iExtIter++;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Take the next time step to tstart + 2*deltat = tstop */
    app->driver->Run(app->iteration_container, app->output, app->integration_container,
                   app->geometry_container, app->solver_container, app->numerics_container,
                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
                   app->interpolator_container, app->transfer_container);

    /* Grab the state vector from su2 */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Get the new Solution from SU2 */
        su2double *appSolution_time_n  = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n();
        su2double *appSolution_time_n1 = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1();
        /* Set the Solution to the xbraid vector */
        u->node[iPoint]->SetSolution(appSolution_time_n);
        u->node[iPoint]->Set_Solution_time_n(appSolution_time_n);
        u->node[iPoint]->Set_Solution_time_n1(appSolution_time_n1);
    }

    /* Grab information about the convergence of the inner iteration */
//    su2double su2Res_rms = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0);
//    if (su2rank == MASTER_NODE)
//    {
//      cout << "rank_t " << braidrank << " moved " << su2size << " processors from " << tstart << " to " << tstop
//           << " with Res_RMS " << su2Res_rms << endl;
//    }

    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    su2double Density_Inf   = app->config_container[ZONE_0]->GetDensity_FreeStreamND();
    su2double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
    su2double Energy_Inf    = app->config_container[ZONE_0]->GetEnergy_FreeStreamND();
    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    CConfig *config         = app->config_container[ZONE_0];

    /* Allocate memory */
    my_Vector* u;
    u          = new my_Vector;
    u->node    = new CVariable*[nPoint];

    /* Initialize the solution vector with the free-stream state */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      u->node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
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

    /* Copy solution from u to v at every Point */
    my_Vector* v;
    v          = new my_Vector;
    v->node    = new CVariable*[nPoint];
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Create new CNSVariable at every Point and initialize with the Solution in u */
        v->node[iPoint]      = new CNSVariable(u->node[iPoint]->GetSolution(),nDim,nVar,config);
        /* Copy Solution at current time n and previous time n-1 from u to v */
        v->node[iPoint]->Set_Solution_time_n(u->node[iPoint]->GetSolution_time_n());
        v->node[iPoint]->Set_Solution_time_n1(u->node[iPoint]->GetSolution_time_n1());
    }

    /* Set the pointer */
    *v_ptr = v;

    return 0;
}

int my_Free( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    /* Delete the CNSVariable at each point in the grid */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        delete u->node[iPoint];
    }

    /* Delete the solution list */
    delete [] u->node;

    /* Delete braid vector */
    delete u;

    return 0;
}

int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Allocate memory for the sum */
    su2double* vec_sum_n  = new su2double[nVar];
    su2double* vec_sum_n1 = new su2double[nVar];

    /* Loop over all points */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Loop over all variables */
        for (int iVar = 0; iVar < nVar; iVar++){
            /* Compute the sum y = alpha x + beta y at time n */
            su2double alphax = alpha * x->node[iPoint]->GetSolution_time_n()[iVar];
            su2double betay  = beta  * y->node[iPoint]->GetSolution_time_n()[iVar];
            vec_sum_n[iVar] = alphax + betay;
            /* Compute the sum y = alpha x + beta y at time n-1 */
            alphax = alpha * x->node[iPoint]->GetSolution_time_n1()[iVar];
            betay  = beta  * y->node[iPoint]->GetSolution_time_n1()[iVar];
            vec_sum_n1[iVar] = alphax + betay;
        }
        /* Store the vector sum in y */
        y->node[iPoint]->Set_Solution_time_n(vec_sum_n);
        y->node[iPoint]->Set_Solution_time_n1(vec_sum_n1);
    }

    /* Destroy vec_sum */
    delete [] vec_sum_n;
    delete [] vec_sum_n1;

    return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute l2norm of the solution list at time n and n1 */
    su2double norm = 0.0;
    /* Loop over all points */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Loop over all variables */
        for (int iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->node[iPoint]->GetSolution_time_n()[iVar], 2);
            norm += pow(u->node[iPoint]->GetSolution_time_n1()[iVar], 2);
        }
    }

    /* Set the pointer */
    *norm_ptr = sqrt(norm);

    return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

    /* Get rank of global communicator */
    int su2rank, braidrank;
    MPI_Comm_rank(app->comm_x, &su2rank);
    MPI_Comm_rank(app->comm_t, &braidrank);

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    /* Retrieve xBraid time information from status object */
    su2double t;
    braid_AccessStatusGetT(astatus, &t);

    /* --- Write Solution_time_n to restart file ---*/

    /* Trick SU2 with the correct iExtIter = (t - t0)/dt - 1  which is used for naming the restart file */
    int iExtIter = (int) round( ( t - app->initialstart ) / app->initialDT) - 1;
//    iExtIter = iExtIter + 10000*braidrank;   /* Make the identifier unique for each braid processor */
    app->config_container[ZONE_0]->SetExtIter(iExtIter);
    /* Check if xBraid tries to write to negative iExtIter */
    if (iExtIter < 0) {
      if (su2rank==MASTER_NODE) cout << "rank_t " << braidrank << " tries to write to iExtIter -1 -> Early Exit.\n" << endl;
      return 0;
    }

    /* Trick SU2 with the current solution for output (SU2 writes CVariable::Solution, not _time_n!) */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        su2double *uSolution = u->node[iPoint]->GetSolution_time_n();
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(uSolution);
    }

    /* Compute the primitive Variables from the conservative ones */
    app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

    /* Call the SU2 output routine */
    if (su2rank==MASTER_NODE) cout << "rank_t " << braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
    app->output->SetResult_Files(app->solver_container, app->geometry_container,
                                 app->config_container, iExtIter, 1);

    /* --- Write Solution_time_n1 to restart file ---*/

    /* Trick SU2 with the correct iExtIter = iExtIter - 1 */
    iExtIter--;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);
    /* Check if xBraid tries to write to negative iExtIter */
    if (iExtIter < 0) {
      if (su2rank==MASTER_NODE) cout << "rank_t " << braidrank << " tries to write to iExtIter -1 -> Early Exit.\n" << endl;
      return 0;
    }

    /* Trick SU2 with the current solution for output */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        su2double *uSolution = u->node[iPoint]->GetSolution_time_n1();
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(uSolution);
    }

    /* Compute the primitive Variables from the conservative ones */
    app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

    /* Call the SU2 output routine */
    if (su2rank==MASTER_NODE) cout << "rank_t " << braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
    app->output->SetResult_Files(app->solver_container, app->geometry_container,
                                 app->config_container, iExtIter, 1);


    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute size of buffer */
    *size_ptr = 2.0 * nPoint * nVar * sizeof(su2double);

    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_Int *size_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Pack the buffer with current and previous time */
    su2double *dbuffer = (su2double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Write Solution at current time to the buffer */
          dbuffer[ibuffer] = u->node[iPoint]->GetSolution_time_n()[iVar];
          ibuffer++;
          /* Write Solution at previous time to the buffer */
          dbuffer[ibuffer] = u->node[iPoint]->GetSolution_time_n1()[iVar];
          ibuffer++;
        }
    }

    /* Compute size of buffer */
    *size_ptr = 2.0 * nPoint * nVar * sizeof(su2double);

    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    su2double Density_Inf   = app->config_container[ZONE_0]->GetDensity_FreeStreamND();
    su2double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
    su2double Energy_Inf    = app->config_container[ZONE_0]->GetEnergy_FreeStreamND();
    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    CConfig *config         = app->config_container[ZONE_0];

    /* Allocate memory for the new braid Vector */
    my_Vector* u;
    u          = new my_Vector;
    u->node    = new CVariable*[nPoint];
    su2double* uSolution    = new su2double[nVar];
    su2double* uSolution_n  = new su2double[nVar];
    su2double* uSolution_n1 = new su2double[nVar];

    /* Initialize the braid vector with the free-stream state */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      u->node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    }

    /* Unpack the buffer and write solution to current and previous time */
    su2double *dbuffer = (su2double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Unpack Solution at current time from the buffer */
          uSolution_n[iVar] = dbuffer[ibuffer];
          ibuffer++;
          /* Unpack Solution at previous time from the buffer */
          uSolution_n1[iVar] = dbuffer[ibuffer];
          ibuffer++;
        }
        /* Write current and previous solution to u*/
        u->node[iPoint]->Set_Solution_time_n(uSolution_n);
        u->node[iPoint]->Set_Solution_time_n1(uSolution_n1);
    }

    /* Delete the uSolution list */
    delete [] uSolution_n;
    delete [] uSolution_n1;

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}
