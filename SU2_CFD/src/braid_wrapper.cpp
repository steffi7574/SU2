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


BraidTape_t* braidTape;

void setupTapeData(){
   /* Create the braid tape */
   braidTape = new BraidTape_t();
}


/* Make a copy of a Vector. Used for primal taping. */
braid_Vector deep_copy( braid_App app, braid_Vector u ){

//    /* Grab variables from the app */
//    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    /* Allocate memory for the new copy v */
//    my_Vector* v = new my_Vector;
//    v->Solution  = new TwoStepSolution(nPoint, nVar);

//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        /* Copy the values from u to v */
//        /* TODO: Set the values with the SetSolution(...), Set_solution_time_n(...) and Set_Solution_timen1(...) */
//        for (int iVar = 0; iVar < nVar; iVar++){
//            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
//            v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
//        }
//    }
//    return v;
}


int my_Step( braid_App        app,
             braid_Vector     ustop,
             braid_Vector     fstop,
             braid_Vector     u,
             braid_StepStatus status ){

//    /* Print action output */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_Step\n");
//    }

//    /* Push the input braid_vector the primal tape */
//    braid_Vector ustore = deep_copy(app, u);
//    braidTape->primal.push_back(ustore);


//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
//    int nDim   = app->geometry_container[ZONE_0][MESH_0]->GetnDim();

//    /* Grab status of current time step from xBraid */
//    double tstart;
//    double tstop;
//    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
////    braid_PhiStatusGetTstartTstop(status, &tstart, &tstop);

//    /* Declare and allocate intermediate casting variables */
//    su2double* cast_n  = new su2double[nVar];
//    su2double* cast_n1 = new su2double[nVar];


//    /* Trick SU2 with xBraid's DeltaT / 2 */
//    double deltat = ( tstop - tstart ) / 2.0;
//    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );

//    /* Trick the su2 solver with the correct state vector (Solution, Solution_time_n and Solution_time_n1*/
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//      for (int iVar = 0; iVar < nVar; iVar++){
//        cast_n[iVar]  = u->Solution->time_n[iPoint][iVar];
//        cast_n1[iVar] = u->Solution->time_n1[iPoint][iVar];
//      }
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast_n);
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(cast_n);
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(cast_n1);
//    }

//    /* Trick SU2 with the correct iExtIter = (t - t0)/dt  -1 */
//    int iExtIter = (int) round( ( tstart + deltat - app->initialstart ) / app->initialDT) ;
//    app->config_container[ZONE_0]->SetExtIter(iExtIter);

//    cout << format("Take 2-step at %d\n", iExtIter);

//    /* Print information output */
//    // if (app->su2rank == MASTER_NODE) cout<<format(" %d: two %1.3f-steps from %1.4f to %1.4f\n", app->braidrank, deltat, tstart, tstop);

//    /* Take the first time step to tstart + deltat */
//    app->driver->Run(app->iteration_container, app->output, app->integration_container,
//                   app->geometry_container, app->solver_container, app->numerics_container,
//                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//                   app->interpolator_container, app->transfer_container);

//    /* Get the objective function */
//    su2double Obj_Func = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag();
//    cout<< format("primal: Obj_Func n1 %1.14e\n", SU2_TYPE::GetValue(Obj_Func));

//    // cout<< format("primal: Obj_Func first step %1.14e\n", SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()));

//    /* Grab the history values from SU2's master node. */
//    if (app->su2rank == MASTER_NODE){
////      /* Grab the flow coefficient values for that time step and store it at time n-1 */
//      u->Solution->Total_CDrag_n1       = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag());
//      u->Solution->Total_CLift_n1       = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift());
//      u->Solution->Total_CSideForce_n1  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CSideForce());
//      u->Solution->Total_CEff_n1        = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CEff());
//      u->Solution->Total_CMx_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMx());
//      u->Solution->Total_CMy_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMy());
//      u->Solution->Total_CMz_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMz());
//      u->Solution->Total_CFx_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFx());
//      u->Solution->Total_CFy_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFy());
//      u->Solution->Total_CFz_n1         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFz());
//    }

//    /*  Store the primal output Braid Vector on the primal tape */
//    my_Vector* ustore_tmp = new my_Vector;
//    ustore_tmp->Solution = new TwoStepSolution(nPoint, nVar);
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//            ustore_tmp->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
//            ustore_tmp->Solution->time_n1[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
//        }
//    }
//    braidTape->primal.push_back(ustore_tmp);

//  //   /* Trick SU2 with the next iExtIter */
//  //   iExtIter++;
//  //   app->config_container[ZONE_0]->SetExtIter(iExtIter);
//  //
//  //   /* Take the next time step to tstart + 2*deltat = tstop */
//  //   // app->driver->Run(app->iteration_container, app->output, app->integration_container,
//  //                 //  app->geometry_container, app->solver_container, app->numerics_container,
//  //                 //  app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//  //                 //  app->interpolator_container, app->transfer_container);
//  //
//  // /* test upddate*/
//  // Obj_Func = 0.0;
//  // for (int iPoint = 0; iPoint < nPoint; iPoint++){
//  //   for (int iDim = 0; iDim < nDim; iDim++){
//  //     /* Compute at Solution */
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(iDim, app->geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord()[iDim] * app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution()[iDim]);
//  //     /* Shift in time */
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n());
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution());
//  //   }
//  //   Obj_Func += app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution()[0];
//  // }
//  // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetTotal_CDrag(Obj_Func);
//  // cout<< format("primal: Obj_Func n1 %1.14e\n", SU2_TYPE::GetValue(Obj_Func));
//  //
//  //   // cout<< format("primal: Obj_Func second step %1.14e\n", SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()));
//  //
//  //   /* Grab the history values from SU2's master node. */
//  //   if (app->su2rank == MASTER_NODE){
//  //     /* Grab the flow coefficient values for that time step and store it at time n-1 */
//  //     u->Solution->Total_CLift_n       = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift());
//  //     u->Solution->Total_CDrag_n       = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag());
//  //     u->Solution->Total_CSideForce_n  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CSideForce());
//  //     u->Solution->Total_CEff_n        = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CEff());
//  //     u->Solution->Total_CMx_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMx());
//  //     u->Solution->Total_CMy_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMy());
//  //     u->Solution->Total_CMz_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CMz());
//  //     u->Solution->Total_CFx_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFx());
//  //     u->Solution->Total_CFy_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFy());
//  //     u->Solution->Total_CFz_n         = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CFz());
//  //   }
//  //

//     /* Grab the solution vectors from su2 for both time steps */
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//            u->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
//            u->Solution->time_n1[iPoint][iVar] = SU2_TYPE::GetValue(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
//        }
//    }

//    /* Free memory of the intermediate casting variables */
//    delete [] cast_n;
//    delete [] cast_n1;

//    /* Tell XBraid no refinement */
//    braid_StepStatusSetRFactor(status, 1);


//    /* Push the braid Action to the action tape */
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::PHI;
//    action->time      = tstart;
//    action->deltat    = deltat;
//    action->StepStatus = status;
//    braidTape->action.push_back(*action);
//    delete action;

//    /*  Store the primal output Braid Vector on the primal tape */
//    braid_Vector ustore_out = deep_copy(app, u);
//    braidTape->primal.push_back(ustore_out);

//    /* Push (->COPY) a pointer to the adjoint Solution list onto the adjoint tape */
//    braidTape->adjoint.push_back(u->Solution_b);


//    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

//    /* Grab variables from the app */
//    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nDim                = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
//    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
//    double Density_Inf   = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetDensity_FreeStreamND());
//    su2double *Velocity_Inf = app->config_container[ZONE_0]->GetVelocity_FreeStreamND();
//    double Energy_Inf    = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetEnergy_FreeStreamND());
//    double Pressure_Inf  = SU2_TYPE::GetValue(app->config_container[ZONE_0]->GetPressure_FreeStreamND());
//    bool compressible = (app->config_container[ZONE_0]->GetKind_Regime() == COMPRESSIBLE);
//    bool incompressible = (app->config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
//    bool freesurface = (app->config_container[ZONE_0]->GetKind_Regime() == FREESURFACE);

//    /* Print action output */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_Init at time %1.4f\n", t);
//    }

//    /* Allocate memory for the primal braid vector */
//    my_Vector* u = new my_Vector;
//    u->Solution  = new TwoStepSolution(nPoint, nVar);

//    /* Set the initial values */
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){

//        /* Initialize the solution with the freestream values */
//        if (compressible) {
//			u->Solution->time_n[iPoint][0]  = Density_Inf;
//			u->Solution->time_n1[iPoint][0] = Density_Inf;
//			for (int iDim = 0; iDim < nDim; iDim++) {
//				u->Solution->time_n[iPoint][iDim+1]  = Density_Inf*SU2_TYPE::GetValue(Velocity_Inf[iDim]);
//				u->Solution->time_n1[iPoint][iDim+1] = Density_Inf*SU2_TYPE::GetValue(Velocity_Inf[iDim]);
//			}
//			u->Solution->time_n[iPoint][nVar-1]  = Density_Inf*Energy_Inf;
//            u->Solution->time_n1[iPoint][nVar-1] = Density_Inf*Energy_Inf;
//		}
//        if (incompressible || freesurface) {
//			u->Solution->time_n[iPoint][0]  = Pressure_Inf;
//			u->Solution->time_n1[iPoint][0] = Pressure_Inf;
//			for (int iDim = 0; iDim < nDim; iDim++) {
//				u->Solution->time_n[iPoint][iDim+1]  = SU2_TYPE::GetValue(Velocity_Inf[iDim])*Density_Inf;
//				u->Solution->time_n1[iPoint][iDim+1] = SU2_TYPE::GetValue(Velocity_Inf[iDim])*Density_Inf;
//			}
//		}
//	}

//    /* Set the pointer */
//    *u_ptr = u;

//    /* Allocate memory for the adjoints of the braid input variables and initialize with zero */
//    (u->Solution_b).reset(new TwoStepSolution(nPoint, nVar));
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//      for (int iVar = 0; iVar < nVar; iVar++){
//        u->Solution_b->time_n[iPoint][iVar] = 0.0;
//      }
//    }
//    cout<< format("u_b %1.14e\n", u->Solution_b->time_n[1][1]);
//    /* Store a pointer to the adjoint on the braid_input tape */
//    braidTape->braid_input_b.push_back(u->Solution_b);


//    /* Push the braid action to the action tape */
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::INIT;
//    braidTape->action.push_back(*action);
//    delete action;

//    return 0;
}

int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr ){

//    /* Grab variables from the app */
//    int nPoint      = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nDim        = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
//    int nVar        = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
//    CConfig *config = app->config_container[ZONE_0];

//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      /* Print information output */
//      if (app->su2rank == MASTER_NODE) cout << format("My_Clone\n");
//    }

//    /* Allocate memory for the new copy v */
//    my_Vector* v = new my_Vector;
//    v->Solution  = new TwoStepSolution(nPoint, nVar);

//    /* Copy the values from u to v */
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
//            v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
//        }
//    }

//    /* Allocate memory for the adjoint SolutionVars and initialize with zero */
//    (v->Solution_b).reset(new TwoStepSolution(nPoint,nVar));
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//      for (int iVar = 0; iVar < nVar; iVar++){
//        v->Solution_b->time_n[iPoint][iVar] = 0.0;
//      }
//    }

//    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
//    braidTape->adjoint.push_back(u->Solution_b);
//    braidTape->adjoint.push_back(v->Solution_b);

//    /* Push the braid action to the action tape*/
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::CLONE;
//    braidTape->action.push_back(*action);
//    delete action;

//    /* Set the pointer */
//    *v_ptr = v;

//    return 0;
}

int my_Free( braid_App app, braid_Vector u ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();

//    /* Print action information */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_Free\n");
//    }

//    /* Delete the primal solution lists (calls destructor of TwoStepSolution) */
//    delete u->Solution;

//    /* Don't delete the adjoint solution list. It is needed in adjoint computation, access it and delete it later via the shared pointer. */

//    /* Delete the braid_Vector */
//    delete u;

//    /* Push the braid action to the action tape */
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::FREE;
//    braidTape->action.push_back(*action);
//    delete action;

//    return 0;
}

int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      /* Print information output */
//      if (app->su2rank == MASTER_NODE) cout << format("My_Sum\n");
//    }

//    /* Compute the sum y = alpha x + beta y at time n and time n-1 */
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//            y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
//                                              + beta  * y->Solution->time_n[iPoint][iVar];
//            y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
//                                              + beta  * y->Solution->time_n1[iPoint][iVar];
//        }
//    }

//    /* Push the braid action */
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::SUM;
//    action->sum_alpha = alpha;
//    action->sum_beta  = beta;
//    braidTape->action.push_back(*action);
//    delete action;

//    /* Push (->COPY) pointers to the adjoint SolVars onto the adjoint tape */
//    braidTape->adjoint.push_back(x->Solution_b);
//    braidTape->adjoint.push_back(y->Solution_b);


//    return 0;
}

int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      /* Print information output */
//      if (app->su2rank == MASTER_NODE) cout << format("My_SpatialNorm\n");
//    }

//    /* Compute l2norm of the solution list at time n and n1 */
//    double norm = 0.0;
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//            norm += pow(u->Solution->time_n[iPoint][iVar], 2);
//            norm += pow(u->Solution->time_n1[iPoint][iVar], 2);
//        }
//    }

//    /* Set the pointer */
//    *norm_ptr = sqrt(norm);

//    return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    /* Allocate memory for the casting variable */
//    su2double* cast = new su2double[nVar];

//    /* Retrieve xBraid time information from status object */
//    double t;
//    braid_AccessStatusGetT(astatus, &t);

//    /* Compute the time step iExtIter = (t - t0)/dt which is used for naming the restart file */
//    int iExtIter = (int) round( ( t - app->initialstart ) / app->initialDT) ;
//    app->config_container[ZONE_0]->SetExtIter(iExtIter);

//    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
//    braidTape->adjoint.push_back(u->Solution_b);

//    /* Print Action Information */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_Access at iExtIter = %d\n", iExtIter);
//    }

//    /* Only continue if iExtIter > 0 !! Otherwise xbraid tries to write at timestep -1*/
//    if (iExtIter>0){


//      /* --- Write Solution_time_n to restart file ---*/

//      /* Trick SU2 with the current solution for output (SU2 writes CVariable::Solution, not _time_n!) */
//      for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0.0; iVar < nVar; iVar++){
//          cast[iVar] = u->Solution->time_n[iPoint][iVar];
//        }
//        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
//      }

//      /* Compute the primitive Variables from the conservative ones */
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

////    /* Call the SU2 output routine */
//      // if (app->su2rank==MASTER_NODE) cout << "rank_t " << app->braidrank << " writes SU2 restart file at iExtIter = " << iExtIter << endl;
//      app->output->SetResult_Files(app->solver_container, app->geometry_container,
//                                 app->config_container, iExtIter, 1);

//      /* Write history values at time n to the app stream */
//      /* NOT WORKING RIGHT NOW....... */
//      //if (app->su2rank == MASTER_NODE){
//      //    *app->history_stream << iExtIter << " " << u->Solution->Total_CLift_n
//      //                         << " " << u->Solution->Total_CDrag_n
//      //                         << " " << u->Solution->Total_CSideForce_n
//      //                         << " " << u->Solution->Total_CMx_n
//      //                         << " " << u->Solution->Total_CMy_n
//      //                         << " " << u->Solution->Total_CMz_n
//      //                         << " " << u->Solution->Total_CFx_n
//      //                         << " " << u->Solution->Total_CFy_n
//      //                         << " " << u->Solution->Total_CFz_n
//      //                         << " " << u->Solution->Total_CEff_n;
//      //    //      for (int iVar = 0; iVar < nVar; iVar++){
//      //    //        *app->history_stream << " " << u->residual_flow_n[iVar];
//      //    //      }
//      //    *app->history_stream << " " << u->residual_dens_n;
//      //    *app->history_stream << "\n";
//      //}


//      /* --- Write Solution_time_n1 to restart file ---*/


//      /* Trick SU2 with the correct iExtIter = iExtIter - 1 */
//      iExtIter--;
//      app->config_container[ZONE_0]->SetExtIter(iExtIter);

//      /* Trick SU2 with the current solution for output */
//      for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0.0; iVar < nVar; iVar++){
//          cast[iVar] = u->Solution->time_n1[iPoint][iVar];
//        }
//        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
//      }

//      /* Compute the primitive Variables from the conservative ones */
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

//      /* Call the SU2 output routine */
//      app->output->SetResult_Files(app->solver_container, app->geometry_container,
//                                       app->config_container, iExtIter, 1);


//      /* Write history values at time n-1 to the app stream */
//      /* NOT WORKING RIGHT NOW....... */
//      //if (app->su2rank == MASTER_NODE){
//      //    *app->history_stream << iExtIter << " " << u->Solution->Total_CLift_n1
//      //                         << " " << u->Solution->Total_CDrag_n1
//      //                         << " " << u->Solution->Total_CSideForce_n1
//      //                         << " " << u->Solution->Total_CMx_n1
//      //                         << " " << u->Solution->Total_CMy_n1
//      //                         << " " << u->Solution->Total_CMz_n1
//      //                         << " " << u->Solution->Total_CFx_n1
//      //                         << " " << u->Solution->Total_CFy_n1
//      //                         << " " << u->Solution->Total_CFz_n1
//      //                         << " " << u->Solution->Total_CEff_n1;
//      //    //      for (int iVar = 0; iVar < nVar; iVar++){
//      //    //        *app->history_stream << " " << u->residual_flow_n1[iVar];
//      //    //      }
//      //    *app->history_stream << " " << u->residual_dens_n1;
//      //    *app->history_stream << "\n";
//      //}

//      /* Add to objective function. */
//      app->Total_Cd_avg += u->Solution->Total_CDrag_n + u->Solution->Total_CDrag_n1;

//    }

//    /* Free memory for the intermediat casting variable */
//    delete [] cast;

//    /* Push the braid action to the action tape*/
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::ACCESS;
//    action->time = t;
//    action->optimiter = app->optimiter;
//    braidTape->action.push_back(*action);
//    delete action;


//    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus  ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    /* Compute size of buffer */
//    *size_ptr = 2.0 * nPoint * nVar * sizeof(double);

//    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus  ){

//    /* Grab variables from the app */
//    int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    /* Print information output */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_BufPack\n");
//    }

//    /* Pack the buffer with current and previous time */
//    double *dbuffer = (double*)buffer;
//    int ibuffer = 0;
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//          /* Write Solution at current time to the buffer */
//          dbuffer[ibuffer] = u->Solution->time_n[iPoint][iVar];
//          ibuffer++;
//          /* Write Solution at previous time to the buffer */
//          dbuffer[ibuffer] = u->Solution->time_n1[iPoint][iVar];
//          ibuffer++;
//        }
//    }

//    /* Set the size of the buffer */
//    int size = 2.0 * nPoint * nVar * sizeof(double);
//    braid_BufferStatusSetSize( bstatus, size );

//    /* Set up the braid action */
//    BraidAction_t* action = new BraidAction_t();
//    action->braidCall = BraidCall_t::BUFPACK;
//    braidTape->action.push_back(*action);
//    delete action;

//    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
//    braidTape->adjoint.push_back(u->Solution_b);


//    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  ){

//    /* Grab variables from the app */
//    int nPoint              = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//    int nVar                = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//    /* Print information output */
//    if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank == MASTER_NODE) cout << format("My_BufUnpack\n");
//    }

//    /* Get the buffer */
//    double *dbuffer = (double*)buffer;
//    int ibuffer = 0;

//    /* Allocate memory for the new braid Vector */
//    my_Vector* u = new my_Vector;
//    u->Solution  = new TwoStepSolution(nPoint, nVar);

//    /* Unpack the buffer and write solution at current and previous time */
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//        for (int iVar = 0; iVar < nVar; iVar++){
//          u->Solution->time_n[iPoint][iVar] = dbuffer[ibuffer];
//          ibuffer++;
//          u->Solution->time_n1[iPoint][iVar] = dbuffer[ibuffer];
//          ibuffer++;
//        }
//    }

//    /* Allocate memory for the adjoint SolutionVars and initialize with zero */
//    (u->Solution_b).reset(new TwoStepSolution(nPoint, nVar));
//    for (int iPoint = 0; iPoint < nPoint; iPoint++){
//      for (int iVar = 0; iVar < nVar; iVar++){
//        u->Solution_b->time_n[iPoint][iVar] = 0.0;
//      }
//    }

//    /* Push (->COPY) the pointer to the adjoint SolVars onto the adjoint tape */
//    braidTape->adjoint.push_back(u->Solution_b);

//    /* Set the pointer */
//    *u_ptr = u;

//    /* Set up the braid action */
//    BraidAction_t* action  = new BraidAction_t();
//    action->braidCall      = BraidCall_t::BUFUNPACK;
//    braidTape->action.push_back(*action);
//    delete action;

//    return 0;
}


void evalAdjointAction( braid_App app, BraidTape_t* braidTape){

//  /* Evaluate the action tape in reverse order */
//  for (std::vector<BraidAction_t>::reverse_iterator action = braidTape->action.rbegin(); action != braidTape->action.rend(); ++action) {
//      switch ( action->braidCall ) {
//        case BraidCall_t::PHI : {
//          my_Step_adjoint(*action, app);
//          break;
//        }
//        case BraidCall_t::INIT : {
//          /* Do nothing. */
//          break;
//        }
//        case BraidCall_t::CLONE : {
//          my_Clone_adjoint(*action, app);
//          break;
//        }
//        case BraidCall_t::FREE : {
//          /* Do nothing. */
//          break;
//        }
//        case BraidCall_t::SUM : {
//          my_Sum_adjoint(*action, app);
//          break;
//        }
//        case BraidCall_t::ACCESS : {
//          my_Access_adjoint(*action, app);
//          break;
//        }
//        case BraidCall_t::BUFPACK : {
//          my_BufPack_adjoint(*action, app);
//          break;
//        }
//        case BraidCall_t::BUFUNPACK : {
//          my_BufUnPack_adjoint(*action, app);
//          break;
//        }
//    }
//  }
}



void my_Step_adjoint( BraidAction_t &action, braid_App app ){

//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank==MASTER_NODE) std::cout << format("%d: PHI adj\n", app->braidrank);
//  }

//  /* Grab variables from the app */
//  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//  int nDim   = app->geometry_container[ZONE_0][MESH_0]->GetnDim();
//  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//  /* Allocate memory for intermediate casting variables */
//  su2double** cast_n  = new su2double*[nPoint];
//  su2double** cast_n1 = new su2double*[nPoint];
//  for(int iPoint=0; iPoint < nPoint; iPoint++){
//    cast_n[iPoint]  = new su2double[nVar];
//    cast_n1[iPoint] = new su2double[nVar];
//  }

//  /* Load the primal vector that was the output  and the input of the xBraid Step function */
//  braid_Vector u_out = braidTape->primal.back();
//  braidTape->primal.pop_back(); /* Pop the pointer */
//  braid_Vector u_tmp = braidTape->primal.back();
//  braidTape->primal.pop_back(); /* Pop the pointer */
//  braid_Vector u_in = braidTape->primal.back();
//  braidTape->primal.pop_back(); /* Pop the pointer */

//  /* Pop the adjoint shared pointer from the tape */
//  std::shared_ptr<TwoStepSolution> usol_b = braidTape->adjoint.back();


//  /* --- Do the same SECOND step that was dont in the in the primal run and record. --- */

//  /* Set the deltat and iExtIter that were used in the primal xbraid run */
//  app->config_container[ZONE_0]->SetDelta_UnstTimeND( action.deltat );
//  int iExtIter = (int) round( (action.time + action.deltat - app->initialstart) / app->initialDT) ;
//  app->config_container[ZONE_0]->SetExtIter(iExtIter+1);

//  // /* Cast the state vector that was used in primal run to su2double and give it to SU2 */
//  // for (int iPoint=0; iPoint < nPoint; iPoint++){
//  //   for (int iVar = 0; iVar < nVar; iVar++){
//  //     cast_n[iVar]  = u_tmp->Solution->time_n[iPoint][iVar];
//  //     cast_n1[iVar] = u_tmp->Solution->time_n1[iPoint][iVar];
//  //   }
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast_n);
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(cast_n);
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(cast_n1);
//  // }
//  //
//  // /* Start CoDi taping. */
//  // AD::StartRecording();
//  //
//  // /* Register input variables */
//  // for (int iPoint=0; iPoint < nPoint; iPoint++){
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->RegisterSolution(true);
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->RegisterSolution_time_n(true);
//  //   app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->RegisterSolution_time_n1(true);
//  // }
//  // app->geometry_container[ZONE_0][MESH_0]->RegisterCoordinates(app->config_container[ZONE_0]);
//  //
//  // /* Tape updating the geometry */
//  // app->geometry_container[ZONE_0][MESH_0]->UpdateGeometry(app->geometry_container[ZONE_0], app->config_container[ZONE_0]);
//  //
//  // /* Record the time step that has been done in primal run */
//  // // app->driver->Run(app->iteration_container, app->output, app->integration_container,
//  //               //  app->geometry_container, app->solver_container, app->numerics_container,
//  //               //  app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//  //               //  app->interpolator_container, app->transfer_container);
//  // /* Get objective function */
//  // // su2double Obj_Func = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag();
//  //
//  // /* test upddate*/
//  // su2double Obj_Func = 0.0;
//  // for (int iPoint = 0; iPoint < nPoint; iPoint++){
//  //   for (int iDim = 0; iDim < nDim; iDim++){
//  //     /* Compute at Solution */
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(iDim, app->geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord()[iDim] * app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution()[iDim]);
//  //     /* Shift in time */
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n());
//  //     app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution());
//  //   }
//  //   Obj_Func += app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution()[0];
//  // }
//  // cout<< format("adjoint: Obj_Func n1 %1.14e\n", SU2_TYPE::GetValue(Obj_Func));
//  //
//  // /* Register Output variables */
//  // if (app->su2rank == MASTER_NODE) AD::RegisterOutput(Obj_Func);
//  //
//  // /* Stop CoDi Taping */
//  // AD::StopRecording();
//  //
//  // /* Initialize the adjoint of the objective function */
//  // if (app->su2rank == MASTER_NODE) SU2_TYPE::SetDerivative(Obj_Func, usol_b->Total_CDrag_n);
//  // usol_b->Total_CDrag_n = 0.0;
//  //
//  // /* Set the adjoint values of the output variables */
//  // for (int iPoint=0; iPoint < nPoint; iPoint++){
//    // for (int iVar = 0; iVar < nVar; iVar++){
//      // cast_n[iVar]  = usol_b->time_n[iPoint][iVar];
//      // cast_n1[iVar] = usol_b->time_n1[iPoint][iVar];
//    // }
//    // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetAdjointSolution_time_n(cast_n);
//    // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetAdjointSolution_time_n1(cast_n1);
//  // }

//  //
//  // /* Evaluate the tape */
//  // AD::ComputeAdjoint();
//  //
//  // /* Get the adjoints from the input variables and store them in the temporary adjoint variable. */
//  // for (int iPoint=0; iPoint < nPoint; iPoint++){
//    /* Use the same order as when registering the input */
//    // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetAdjointSolution(app->tmpadj[iPoint]);
//    // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetAdjointSolution_time_n(app->tmpadj_n[iPoint]);
//    // app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetAdjointSolution_time_n1(app->tmpadj_n1[iPoint]);
//  // }
//  //
//  // /* Get the reduced gradient. */
//  // su2double *Coord;
//  // for (int iPoint = 0; iPoint < nPoint; iPoint++){
//  //   Coord = app->geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord();
//  //   for (int iDim=0; iDim < nDim; iDim++){
//  //     double sens = SU2_TYPE::GetDerivative(Coord[iDim]);
//  //     if (iPoint == 8737 && iDim == 0) cout<< format("sens %1.14e\n", sens);
//  //     app->redgrad[iPoint][iDim] += sens;
//  //   }
//  // }
//  //
//  //
//  // /* Reset the CoDi Tape */
//  // AD::Reset();
//  //



//  /* --- Do the same FIRST step that was done in the in the primal run and record. --- */

//  /* Reset the CoDi Tape */
//  AD::Reset();

//  /* Set the deltat and iExtIter that were used in the primal xbraid run */
//  app->config_container[ZONE_0]->SetExtIter(iExtIter);

//  /* Cast the state vector that was used in primal run to su2double */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      cast_n[iPoint][iVar]  = u_in->Solution->time_n[iPoint][iVar];
//      cast_n1[iPoint][iVar] = u_in->Solution->time_n1[iPoint][iVar];
//    }
//  }

//  /* Start CoDi taping. */
//  AD::StartRecording();

//  /* Register input variables */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      AD::globalTape.registerInput(cast_n[iPoint][iVar]);
//      AD::globalTape.registerInput(cast_n1[iPoint][iVar]);
//    }
//  }
//  /* Register coordinates */
//  for (int iPoint = 0; iPoint < nPoint; iPoint++){
//      for (int iDim = 0; iDim < nDim; iDim ++){
//        AD::globalTape.registerInput(app->geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord()[iDim]);
//    }
//  }

//  /* Record updating the geometry */
//  app->geometry_container[ZONE_0][MESH_0]->UpdateGeometry(app->geometry_container[ZONE_0], app->config_container[ZONE_0]);

//  /* Set the solution from the casting variables */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast_n[iPoint]);
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(cast_n[iPoint]);
//      app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(cast_n1[iPoint]);
//  }


//  /* Record an update */
//  app->driver->Run(app->iteration_container, app->output, app->integration_container,
//                   app->geometry_container, app->solver_container, app->numerics_container,
//                   app->config_container, app->surface_movement, app->grid_movement, app->FFDBox,
//                   app->interpolator_container, app->transfer_container);

//  /* Get the objective function */
//  su2double Obj_Func = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag();
//  cout<< format("adjoint: Obj_Func n1 %1.14e\n", SU2_TYPE::GetValue(Obj_Func));

//  /* Register Output variables */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      AD::globalTape.registerOutput(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
//      AD::globalTape.registerOutput(app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
//    }
//  }
//  AD::globalTape.registerOutput(Obj_Func);

//  /* Stop CoDi Taping */
//  AD::StopRecording();

//  /* Initialize the adjoint of the objective function */
//  // Obj_Func.setGradient(1.0);
// Obj_Func.setGradient(usol_b->Total_CDrag_n1);
//  // if (app->su2rank == MASTER_NODE) SU2_TYPE::SetDerivative(Obj_Func, usol_b->Total_CDrag_n1);
//  usol_b->Total_CDrag_n1 = 0.0;

//  /* Set the adjoint values of the output variables to the temporary adjoints */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      (app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]).setGradient(usol_b->time_n[iPoint][iVar]);
//      (app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]).setGradient(usol_b->time_n1[iPoint][iVar]);
//    }
//  }
//  /* Evaluate the tape */
//  AD::ComputeAdjoint();

//  /* Get the adjoints from the input variables and store them in the xbraid adjoints. */
//  for (int iPoint=0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      usol_b->time_n[iPoint][iVar] = (cast_n[iPoint][iVar]).getGradient();
//      usol_b->time_n1[iPoint][iVar] = (cast_n1[iPoint][iVar]).getGradient();
//    }
//  }

//  /* Get the test reduced gradient. */
//  for (int iPoint = 0; iPoint < nPoint; iPoint++){
//    for (int iDim = 0; iDim < nDim; iDim++){
//      double localsens = (app->geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord()[iDim]).getGradient();
//      if (iPoint == 8737 ) cout<< format("sens %1.14e\n", localsens);
//      app->redgrad[iPoint][iDim] += localsens;
//    }
//  }

//  // cout<< format("Test usol_b out %d %d %1.14e %1.14e\n", iPoint, iVar, usol_b->time_n[iPoint][iVar], usol_b->time_n1[iPoint][iVar]);


//  /* Free memory of the intermediate casting vectors */
//  for (int iPoint = 0; iPoint < nPoint; iPoint++){
//    delete [] cast_n[iPoint];
//    delete [] cast_n1[iPoint];
//  }

//  /* Free the memory of the intermediate primal vectors. */
//  delete u_out->Solution;
//  delete u_out;
//  delete u_tmp->Solution;
//  delete u_tmp;
//  delete u_in->Solution;
//  delete u_in;

//  /* Pop the adjoint pointer from the vector */
//  (braidTape->adjoint).pop_back();
//  usol_b.reset();

}


void my_Access_adjoint( BraidAction_t &action , braid_App app ){
/*
 * my_Access:              Total_Cd_avg += Total_CDrag_n + Total_CDrag_n1
 * my_Access_adjoint:  Total_CDrag_nb   += Total_Cd_avg_b
 *                     Total_CDrag_n1_b += Total_Cd_avg_b

  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
      if (app->su2rank==MASTER_NODE) std::cout << format("%d: Access adj\n", app->braidrank);
  }

//  /* Grab variables from the app */
//  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();
//  /* Compute the time step iExtIter = (t - t0)/dt which is used for naming the restart file */
//  int iExtIter = (int) round( ( action.time - app->initialstart ) / app->initialDT) ;

//  /* Print action information */
//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//    if (app->su2rank == MASTER_NODE) cout << format("My_Access_adj at iExtIter = %d\n", iExtIter);
//  }


//  /* Pop the adjoint from the tape */
//  std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
//  (braidTape->adjoint).pop_back();


//  /* If CPoint: Set the adjoint seed from previous iteration */
//  int cfactor = app->config_container[ZONE_0]->GetBraid_CFactor();
//  if (app->ncpoints > 1){
//    if (_braid_IsCPoint(iExtIter/2,cfactor) ) // Divide by two because 2-Step XBraid !
//    { /* ATTENTION: THIS IS A BUGGY, IF MINCOARSE FORCES TO DO SERIAL RUN, then ncpoint = 1 */

//      /* Store the pointer in the braid_output_b vector */
//      int pos = (int) (iExtIter/2 - app->ilower) / cfactor;
//      braidTape->braid_output_b[pos] = usol_b;

//      /* Set the adjoint seed */
//      for (int iPoint=0; iPoint < nPoint; iPoint++){
//        for (int iVar=0; iVar < nVar; iVar++){
//          usol_b->time_n[iPoint][iVar]  = app->optimadjoint[pos]->time_n[iPoint][iVar];
//          usol_b->time_n1[iPoint][iVar] = app->optimadjoint[pos]->time_n1[iPoint][iVar];
//        }
//      }
//    }
//  }

//  /* Set the seed for the costfunction */
//  if (iExtIter>0){
//    usol_b->Total_CDrag_n  += app->Total_Cd_avg_b;
//    usol_b->Total_CDrag_n1 += app->Total_Cd_avg_b;
//  }

//  /* Delete the shared pointer */
//  usol_b.reset();

}


void my_Sum_adjoint( BraidAction_t &action, braid_App app ){
  /*
  * my_Sum:        y = alpha x + beta y;
  * my_sum_adj:     x_b += alpha y_b;
  *                 y_b = beta y_b;
  */

//  /* Print action information */
//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank==MASTER_NODE) std::cout << format("%d: SUM adj\n", app->braidrank);
//  }

//  /* Grab variables from the app */
//  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//  /* Get the coefficients of the sum from the action */
//  double alpha = action.sum_alpha;
//  double beta  = action.sum_beta;

//  /* Pop the adjoint from the tape */
//  std::shared_ptr<TwoStepSolution> ysol_b = (braidTape->adjoint).back();
//  (braidTape->adjoint).pop_back();
//  std::shared_ptr<TwoStepSolution> xsol_b = (braidTape->adjoint).back();
//  (braidTape->adjoint).pop_back();

//  /* Perform adjoint steps for all points */
//  for (int iPoint = 0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      /* Time n */
//      xsol_b->time_n[iPoint][iVar] += alpha * ysol_b->time_n[iPoint][iVar];
//      ysol_b->time_n[iPoint][iVar] = beta * ysol_b->time_n[iPoint][iVar];
//      /* Time n-1 */
//      xsol_b->time_n1[iPoint][iVar] += alpha * ysol_b->time_n1[iPoint][iVar];
//      ysol_b->time_n1[iPoint][iVar] = beta * ysol_b->time_n1[iPoint][iVar];
//    }
//  }

//  /* Delete the shared pointer */
//  xsol_b.reset();
//  ysol_b.reset();

}


void my_Clone_adjoint( BraidAction_t &action, braid_App app ){
  /* my_Clone:        v  = u;
   * my_Clone_adj:  u_b += v_b;
   *                v_b  = 0
   */

//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank==MASTER_NODE) std::cout << format("%d: Clone adj\n", app->braidrank);
//  }

//  /* Pop the adjoint from the tape */
//  std::shared_ptr<TwoStepSolution> vsol_b = (braidTape->adjoint).back();
//  (braidTape->adjoint).pop_back();
//  std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
//  (braidTape->adjoint).pop_back();

//  /* Grab variables from the app */
//  int nPoint = app->geometry_container[ZONE_0][MESH_0]->GetnPoint();
//  int nVar   = app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar();

//  /* Perform adjoint steps for all points */
//  for (int iPoint = 0; iPoint < nPoint; iPoint++){
//    for (int iVar = 0; iVar < nVar; iVar++){
//      /* Time n */
//      usol_b->time_n[iPoint][iVar] += vsol_b->time_n[iPoint][iVar];
//      vsol_b->time_n[iPoint][iVar] = 0.0;
//      /* Time n-1 */
//      usol_b->time_n1[iPoint][iVar] += vsol_b->time_n1[iPoint][iVar];
//      vsol_b->time_n1[iPoint][iVar] = 0.0;
//    }
//  }

//  /* Delete the shared pointer */
//  usol_b.reset();
//  vsol_b.reset();
}


void my_BufPack_adjoint( BraidAction_t &action, braid_App app ){
  /* my_Bufpack:       buffer = u;
 *                   MPI_Isend(buffer);
 * my_BufPack_adj:   MPI_recv(buffer);
 *                   u_b = buffer;    */
//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank==MASTER_NODE) std::cout << format("%d: BufPack adj\n", app->braidrank);
//  }

// /* Pop the adjoint from the tape */
// std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
// (braidTape->adjoint).pop_back();

// /* TODO: Receive and Unpack the buffer into the adjoint vars */

// /* Remove the shared pointer */
// usol_b.reset();
}

void my_BufUnPack_adjoint( BraidAction_t &action, braid_App app ){
  /* my_BufUnpack:       MPI_Irecv (buffer);
   *                     u = buffer;
   * my_BufUnpack_adj:   buffer = u_b;
   *                     MPI_send(buffer);   */

//  if (app->config_container[ZONE_0]->GetBraid_Action_Verb()){
//      if (app->su2rank==MASTER_NODE) std::cout << format("%d: BufUnPack adj\n", app->braidrank);
//  }

//   /* Pop the adjoint from the tape */
//   std::shared_ptr<TwoStepSolution> usol_b = (braidTape->adjoint).back();
//   (braidTape->adjoint).pop_back();

//   /* Remove the shared pointer */
//   usol_b.reset();


}
