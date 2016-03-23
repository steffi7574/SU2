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

    /* Allocate memory on the heap */
    my_Vector* u;
    u = new my_Vector;

    /* Initialize the vector with initial values
    // from Constructor of CNSSolver(geom, config, imglevel) :
    node = new CVariable*[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
    */

    /* Set pointer */
    *u_ptr = u;


  return 0;
}

int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr ){

  return 0;
}

int my_Free( braid_App app, braid_Vector u ){

  return 0;
}

int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y ){

  return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

  return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

  return 0;
}

int my_BufSize ( braid_App app, int *size_ptr ){

  return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer,
                braid_Int *size_ptr ){

  return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr ){

  return 0;
}
