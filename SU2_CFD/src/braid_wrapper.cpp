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


BraidTape_t* braidTape;

void setupTapeData(){
   /* Create the braid tape */
   braidTape = new BraidTape_t();
}


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

    /* Print action output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_Step\n");
    }

    /* Push the input braid_vector the primal tape */
    braid_Vector ustore = deep_copy(app, u);
    braidTape->primal.push_back(ustore);


    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Grab status of current time step from xBraid */
    double tstart;
    double tstop;
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
//    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

    /* Trick SU2 with xBraid's DeltaT / 2 */
    double deltat = ( tstop - tstart ) / 2.0;
    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );

    /* Trick the su2 solver with the correct state vector (Solution, Solution_time_n and Solution_time_n1*/
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution->time_n[iPoint]);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(u->Solution->time_n[iPoint]);
        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(u->Solution->time_n1[iPoint]);
    }

    /* Trick SU2 with the correct iExtIter = (t - t0)/dt  -1 */
    int iExtIter = (int) round( ( tstart + deltat - app->initialstart ) / app->initialDT) ;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Print information output */
    // if (app->su2rank == MASTER_NODE) cout<<format(" %d: two %1.3f-steps from %1.4f to %1.4f\n", app->braidrank, deltat, tstart, tstop);

    /* Take the first time step to tstart + deltat */
    app->driver->Run(app->iteration_container, app->output, app->integration_container,
                   app->geometry_container, app->solver_container, app->numerics_container,
                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
                   app->interpolator_container, app->transfer_container);

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
    app->driver->Run(app->iteration_container, app->output, app->integration_container,
                   app->geometry_container, app->solver_container, app->numerics_container,
                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
                   app->interpolator_container, app->transfer_container);

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
            u->Solution->time_n[iPoint][iVar]  = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar];
            u->Solution->time_n1[iPoint][iVar] = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar];
        }
    }

    /* Tell XBraid no refinement */
    braid_StepStatusSetRFactor(status, 1);


    /* Push the braid Action to the action tape */
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::PHI;
    action->deltat    = deltat;
    action->StepStatus = status;
    braidTape->action.push_back(*action);
    delete action;

    /*  Store the primal output Braid Vector on the primal tape */
    braid_Vector ustore_out = deep_copy(app, u);
    braidTape->primal.push_back(ustore_out);

    /* Push (->COPY) a pointer to the adjoint Solution list onto the adjoint tape */
    braidTape->adjoint.push_back(u->Solution_b);


    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
    double Density_Inf   = app->config_container[ZONE_0]->GetDensity_FreeStreamND();
    double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
    double Energy_Inf    = app->config_container[ZONE_0]->GetEnergy_FreeStreamND();
    double Pressure_Inf  = app->config_container[ZONE_0]->GetPressure_FreeStreamND();
    bool compressible = (app->config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (app->config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
    bool freesurface = (app->config_container[ZONE_0]->GetKind_Regime() == FREESURFACE);

    /* Print action output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_Init at time %1.4f\n", t);
    }

    /* Allocate memory for the primal braid vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(nPoint, nVar);

    /* Set the initial values */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){

        /* Initialize the solution with the freestream values */
        if (compressible) {
			u->Solution->time_n[iPoint][0]  = Density_Inf;
			u->Solution->time_n1[iPoint][0] = Density_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution->time_n[iPoint][iDim+1]  = Density_Inf*Velocity_Inf[iDim];
				u->Solution->time_n1[iPoint][iDim+1] = Density_Inf*Velocity_Inf[iDim];
			}
			u->Solution->time_n[iPoint][nVar-1]  = Density_Inf*Energy_Inf;
            u->Solution->time_n1[iPoint][nVar-1] = Density_Inf*Energy_Inf;
		}
        if (incompressible || freesurface) {
			u->Solution->time_n[iPoint][0]  = Pressure_Inf;
			u->Solution->time_n1[iPoint][0] = Pressure_Inf;
			for (int iDim = 0; iDim < nDim; iDim++) {
				u->Solution->time_n[iPoint][iDim+1]  = Velocity_Inf[iDim]*Density_Inf;
				u->Solution->time_n1[iPoint][iDim+1] = Velocity_Inf[iDim]*Density_Inf;
			}
		}
	}

    /* Set the pointer */
    *u_ptr = u;

    /* Allocate memory for the adjoints of the braid input variables and initialize with zero */
    (u->Solution_b).reset(new TwoStepSolution(nPoint, nVar));
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      for (int iVar = 0; iVar < nVar; iVar++){
        u->Solution_b->time_n[iPoint][iVar] = 0.0;
      }
    }
    cout<< format("u_b %1.14e\n", u->Solution_b->time_n[1][1]);
    /* Store a pointer to the adjoint on the braid_input tape */
    braidTape->braid_input_b.push_back(u->Solution_b);


    /* Push the braid action to the action tape */
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::INIT;
    braidTape->action.push_back(*action);
    delete action;

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
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(nPoint, nVar);

    /* Copy the values from u to v */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
        }
    }

    /* Allocate memory for the adjoint SolutionVars and initialize with zero */
    (v->Solution_b).reset(new TwoStepSolution(nPoint,nVar));
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      for (int iVar = 0; iVar < nVar; iVar++){
        v->Solution_b->time_n[iPoint][iVar] = 0.0;
      }
    }

    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
    braidTape->adjoint.push_back(u->Solution_b);
    braidTape->adjoint.push_back(v->Solution_b);

    /* Push the braid action to the action tape*/
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::CLONE;
    braidTape->action.push_back(*action);
    delete action;

    /* Set the pointer */
    *v_ptr = v;

    return 0;
}

int my_Free( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

    /* Print action information */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_Free\n");
    }

    /* Delete the primal solution lists (calls destructor of TwoStepSolution) */
    delete u->Solution;

    /* Don't delete the adjoint solution list. It is needed in adjoint computation, access it and delete it later via the shared pointer. */

    /* Delete the braid_Vector */
    delete u;

    /* Push the braid action to the action tape */
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::FREE;
    braidTape->action.push_back(*action);
    delete action;

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

    /* Compute the sum y = alpha x + beta y at time n and time n-1 */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
                                              + beta  * y->Solution->time_n[iPoint][iVar];
            y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
                                              + beta  * y->Solution->time_n1[iPoint][iVar];
        }
    }

    /* Push the braid action */
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::SUM;
    action->sum_alpha = alpha;
    action->sum_beta  = beta;
    braidTape->action.push_back(*action);
    delete action;

    /* Push (->COPY) pointers to the adjoint SolVars onto the adjoint tape */
    braidTape->adjoint.push_back(x->Solution_b);
    braidTape->adjoint.push_back(y->Solution_b);


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

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Retrieve xBraid time information from status object */
    double t;
    braid_AccessStatusGetT(astatus, &t);

    /* Compute the time step iExtIter = (t - t0)/dt which is used for naming the restart file */
    int iExtIter = (int) round( ( t - app->initialstart ) / app->initialDT) ;
    app->config_container[ZONE_0]->SetExtIter(iExtIter);

    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
    braidTape->adjoint.push_back(u->Solution_b);


    /* Only continue if iExtIter > 0 !! Otherwise xbraid tries to write at timestep -1*/
    if (iExtIter>0){

       /* Print Action Information */
       if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
          if (app->su2rank == MASTER_NODE) cout << format("My_Access at iExtIter = %d\n", iExtIter);
       }

      /* Push the input braid_vector the primal tape */
      braid_Vector ustore = deep_copy(app, u);
      braidTape->primal.push_back(ustore);


      /* --- Write Solution_time_n to restart file ---*/

      /* Trick SU2 with the current solution for output (SU2 writes CVariable::Solution, not _time_n!) */
      for (int iPoint = 0; iPoint < nPoint; iPoint++){
          app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution->time_n[iPoint]);
      }

      /* Compute the primitive Variables from the conservative ones */
      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

//    /* Call the SU2 output routine */
      // if (app->su2rank==MASTER_NODE) cout << "rank_t " << app->braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
      app->output->SetResult_Files(app->solver_container, app->geometry_container,
                                 app->config_container, iExtIter, 1);

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

      /* Trick SU2 with the current solution for output */
      for (int iPoint = 0; iPoint < nPoint; iPoint++){
          app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(u->Solution->time_n1[iPoint]);
      }

      /* Compute the primitive Variables from the conservative ones */
      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

      /* Call the SU2 output routine */
      app->output->SetResult_Files(app->solver_container, app->geometry_container,
                                       app->config_container, iExtIter, 1);


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

      /* Add to objective function. */
      app->Total_Cd_avg += u->Total_CDrag_n + u->Total_CDrag_n1;

    }

    /* Push the braid action to the action tape*/
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::ACCESS;
    action->time = t;
    action->optimiter = app->optimiter;
    braidTape->action.push_back(*action);
    delete action;


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

    /* Print information output */
    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank == MASTER_NODE) cout << format("My_BufPack\n");
    }

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

    /* Set up the braid action */
    BraidAction_t* action = new BraidAction_t();
    action->braidCall = BraidCall_t::BUFPACK;
    braidTape->action.push_back(*action);
    delete action;

    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
    braidTape->adjoint.push_back(u->Solution_b);


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

    /* Allocate memory for the adjoint SolutionVars and initialize with zero */
    (u->Solution_b).reset(new TwoStepSolution(nPoint, nVar));
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
      for (int iVar = 0; iVar < nVar; iVar++){
        u->Solution_b->time_n[iPoint][iVar] = 0.0;
      }
    }

    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
    braidTape->adjoint.push_back(u->Solution_b);

    /* Set the pointer */
    *u_ptr = u;

    /* Set up the braid action */
    BraidAction_t* action  = new BraidAction_t();
    action->braidCall      = BraidCall_t::BUFUNPACK;
    braidTape->action.push_back(*action);
    delete action;

    return 0;
}


void evalAdjointAction( braid_App app, BraidTape_t* braidTape){

  /* Evaluate the action tape in reverse order */
  for (std::vector<BraidAction_t>::reverse_iterator action = braidTape->action.rbegin(); action != braidTape->action.rend(); ++action) {
      switch ( action->braidCall ) {
        case BraidCall_t::PHI : {
          my_Step_adjoint(*action, app);
          break;
        }
        case BraidCall_t::INIT : {
          /* Do nothing. */
          break;
        }
        case BraidCall_t::CLONE : {
          my_Clone_adjoint(*action, app);
          break;
        }
        case BraidCall_t::FREE : {
          /* Do nothing. */
          break;
        }
        case BraidCall_t::SUM : {
          my_Sum_adjoint(*action, app);
          break;
        }
        case BraidCall_t::ACCESS : {
          my_Access_adjoint(*action, app);
          break;
        }
        case BraidCall_t::BUFPACK : {
          my_BufPack_adjoint(*action, app);
          break;
        }
        case BraidCall_t::BUFUNPACK : {
          my_BufUnPack_adjoint(*action, app);
          break;
        }
    }
  }
}



void my_Step_adjoint( BraidAction_t &action, braid_App app ){

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: PHI adj\n", app->braidrank);
  }

  /* Grab variables from the app */
  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

  /* Load the primal vector that was the output  and the input of the xBraid Step function */
  braid_Vector u_out = braidTape->primal.back();
  braidTape->primal.pop_back(); /* Pop the pointer */
  braid_Vector u_in = braidTape->primal.back();
  braidTape->primal.pop_back(); /* Pop the pointer */

  /* Pop the adjoint shared pointer from the tape */
  std::shared_ptr<TwoStepSolution> usol_b = braidTape->adjoint.back();

  /* Set the Time step that was used in the primal xbraid run */
  app->config_container[ZONE_0]->SetDelta_UnstTimeND( action.deltat );

  /* TODO: Implement adjoint action */


  /* Free the memory of the intermediate primal vectors. */
  delete u_out->Solution;
  delete u_out;
  delete u_in->Solution;
  delete u_in;

  /* Pop the adjoint pointer from the vector */
  (braidTape->adjoint).pop_back();
  usol_b.reset();

}


void my_Access_adjoint( BraidAction_t &action , braid_App app ){

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: Access adj\n", app->braidrank);
  }

  /* Grab variables from the app */
  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
  /* Compute the time step iExtIter = (t - t0)/dt which is used for naming the restart file */
  int iExtIter = (int) round( ( action.time - app->initialstart ) / app->initialDT) ;

  /* Pop the adjoint from the tape */
  std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
  (braidTape->adjoint).pop_back();


  /* TODO: If CPoint: Set the adjoint seed from previous iteration */

  if (iExtIter>0){

    /* Load the primal braid vector that was used in the primal xbraid run */
    my_Vector* u_in = braidTape->primal.back();
    braidTape->primal.pop_back();
    /* TODO: Implement adjoint of cost function */

    /* Free the memory of the intermediate primal vector. */
    delete u_in->Solution;
    delete u_in;
  }

  /* Delete the shared pointer */
  usol_b.reset();

}


void my_Sum_adjoint( BraidAction_t &action, braid_App app ){

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: SUM adj\n", app->braidrank);
  }

  /* Get the coefficients of the sum */
  double alpha = action.sum_alpha;
  double beta  = action.sum_beta;

  /* Pop the adjoint from the tape */
  std::shared_ptr<TwoStepSolution> ysol_b = (braidTape->adjoint).back();
  (braidTape->adjoint).pop_back();
  std::shared_ptr<TwoStepSolution> xsol_b = (braidTape->adjoint).back();
  (braidTape->adjoint).pop_back();

  /* TODO: Implement the adjoint action  */

  /* Delete the shared pointer */
  xsol_b.reset();
  ysol_b.reset();

}


void my_Clone_adjoint( BraidAction_t &action, braid_App app ){
  /* my_Clone:        v  = u;
   * my_Clone_adj:  u_b += v_b;
   *                v_b  = 0      */

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: Clone adj\n", app->braidrank);
  }

  /* Pop the adjoint from the tape */
  std::shared_ptr<TwoStepSolution> vsol_b = (braidTape->adjoint).back();
  (braidTape->adjoint).pop_back();
  std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
  (braidTape->adjoint).pop_back();


  /* TODO: Implement the adjoint actions */

  /* Remove the shared pointer */
  usol_b.reset();
  vsol_b.reset();
}


void my_BufPack_adjoint( BraidAction_t &action, braid_App app ){
  /* my_Bufpack:       buffer = u;
 *                   MPI_Isend(buffer);
 * my_BufPack_adj:   MPI_recv(buffer);
 *                   u_b = buffer;    */
  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: BufPack adj\n", app->braidrank);
  }

 /* Pop the adjoint from the tape */
 std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
 (braidTape->adjoint).pop_back();

 /* TODO: Receive and Unpack the buffer into the adjoint vars */

 /* Remove the shared pointer */
 usol_b.reset();
}

void my_BufUnPack_adjoint( BraidAction_t &action, braid_App app ){
  /* my_BufUnpack:       MPI_Irecv (buffer);
   *                     u = buffer;
   * my_BufUnpack_adj:   buffer = u_b;
   *                     MPI_send(buffer);   */

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: BufUnPack adj\n", app->braidrank);
  }

   /* Pop the adjoint from the tape */
   std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
   (braidTape->adjoint).pop_back();

   /* Remove the shared pointer */
   usol_b.reset();


}
