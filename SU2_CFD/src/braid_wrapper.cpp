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

    /* Only implemented for single Zone problems */
    int iZone = 0;


    /* Grab status of current time step */
    su2double tstart;
    su2double tstop;
    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

    /* Trick the su2 solver with the new DeltaT */
    app->config_container[iZone]->SetDelta_UnstTimeND(tstop-tstart);

    /* Trick the su2 solver with the right state vector */
//    app->solver_containter[iZone][iMGLevel]->node = u->node;
//    ODER besser:
//      forall iPoints:
//        app->solver_containter[iZone][iMGLevel]->node[iPoint]->SetSolution(u->node->GetSolution())


    /* Take a time step */
    // app->driver->Run(app->iteration_container, app->output, app->integration_container,
                  // app->geometry_container, app->solver_container, app->numerics_container,
                  // app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
                  // app->interpolator_container, app->transfer_container);

    /* Grab the state vector from su2 */
//    forall iPoints
//        u->node[iPoint]->SetSolution(app->solver_container[iZone][iMGLevel]->node[iPoint]->GetSolution())

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
        su2double *uSolution = u->node[iPoint]->GetSolution();
        v->node[iPoint]      = new CNSVariable(uSolution,nDim,nVar,config);
        /* Copy Solution at previous time n from u to v */
        v->node[iPoint]->SetSolution_time_n(u->node[iPoint]->GetSolution_time_n());
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

    /* Loop over all points */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Loop over all variables */
        for (int iVar = 0; iVar < nVar; iVar++){
            /* Compute the sum y = alpha * x + beta * y  */
            su2double alphax = alpha * x->node[iPoint]->GetSolution(iVar);
            su2double betay  = beta  * y->node[iPoint]->GetSolution(iVar);
            y->node[iPoint]->SetSolution(iVar, alphax + betay);
        }
    }

    return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute l2norm of the solution list */
    su2double norm = 0.0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->node[iPoint]->GetSolution(iVar), 2);
        }
    }

  /* Set the pointer */
  *norm_ptr = sqrt(norm);

  return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nDim   = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Retrieve xBraid time information from status object */
    su2double t;
    braid_AccessStatusGetT(astatus, &t);

    /* Compute the current iExtIter for naming the output file and pass it to SU2 */
    int iExtIter = ( t - app->tstart ) / app->initialDT ;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Trick SU2 with the current state / solution */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        su2double *uSolution = u->node[iPoint]->GetSolution();
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(uSolution);
    }

    /* Call the SU2 output routine */
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
          dbuffer[ibuffer] = u->node[iPoint]->GetSolution(iVar);
          ibuffer++;
          /* Write Solution at previous time to the buffer */
          dbuffer[ibuffer] = u->node[iPoint]->GetSolution_time_n(iVar);
          ibuffer++;
        }
    }

    /* Compute size of buffer */
    *size_ptr = 2.0 * nPoint * nVar * sizeof(su2double);

    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr ){

    /* Initialize new braid vector */
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


    /* Unpack the buffer and write solution to current and previous time */
    su2double *dbuffer = (su2double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
          /* Write Solution at current time from the buffer */
          u->node[iPoint]->SetSolution(iVar,dbuffer[ibuffer]);
          ibuffer++;
          /* Write Solution at previous time from the buffer */
          u->node[iPoint]->SetSolution_time_n(iVar,dbuffer[ibuffer]);
          ibuffer++;
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}
