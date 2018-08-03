/*!
 * \file braid_wrapper.cpp
 * \brief Functions for XBraid integration.
 *
 * \author S. Guenther
 *
 */


//#include <braid.hpp>
#include "../include/util.hpp"
#include "../include/braid_structure.hpp"
#include "../../Common/include/ad_structure.hpp"


/* Make a copy of a Vector. Used for primal taping. */
braid_Vector deep_copy( braid_App app, braid_Vector u ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();

    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        /* Copy the values from u to v */
        /* TODO: Set the values with the SetSolution(...), Set_solution_time_n(...) and Set_Solution_timen1(...) */
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            if(app->BDF2) v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
        }
    }
    return v;
}


int my_Step( braid_App        app,
             braid_Vector     ustop,
             braid_Vector     fstop,
             braid_Vector     u,
             braid_StepStatus status){

    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    
    double  tstart, tstop, deltat;
    unsigned long nPoint;
    unsigned long  iExtIter;
    double  abs_accuracy, rel_accuracy;
    unsigned short  FinestMesh = app->config->GetFinestMesh(),  nVar, nVar_Turb = 0, iVar;
    int level;
    /* Set the SU2 accuracy */
    app->config->SetOrderMagResidual(app->SU2_OrderMagResidual);
    app->config->SetMinLogResidual(app->SU2_MinLogResidual);

    /* Reduce accuracy on coarser grid levels */
    braid_StepStatusGetLevel(status, &level);
    if (level > 0 ){
        rel_accuracy = app->config->GetBraid_CoarsegridAccur_rel();
        if (rel_accuracy > 0) app->config->SetOrderMagResidual(rel_accuracy);
        abs_accuracy = app->config->GetBraid_CoarsegridAccur_abs();
        if (abs_accuracy < 0) app->config->SetMinLogResidual(abs_accuracy);
        
        if (app->config->GetBraid_CoarseGrid_Space()){ 
          app->config->SetFinestMesh(MESH_0+1);
        }
    }

    /* Grab variables from the app */
    nPoint = app->geometry[MESH_0]->GetnPoint();
    nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }

    /* Set the time-step size and time-step index */
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
    if (app->BDF2) {
        deltat = (tstop - tstart) / 2.0; 
    } else {
        deltat = (tstop - tstart);
    }
    iExtIter = (int) round( tstop  / app->initialDT );

    app->config->SetDelta_UnstTimeND( deltat );
    app->config->SetExtIter(iExtIter); //TODO: Check if BDF2

    /* Trick SU2 with the correct solution values */ 
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        app->solver[MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(iVar, u->Solution->time_n[iPoint][iVar]);
      }
      if (turbulent){
        for (iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
          app->solver[MESH_0][TURB_SOL]->node[iPoint]->SetSolution(iVar-nVar, u->Solution->time_n[iPoint][iVar]);
        }    
      }
    }
    

    if (turbulent){
     
      /* Compute primitive variables (needed for turbulent post-processing) */
      
      app->solver[MESH_0][FLOW_SOL]->Preprocessing(app->geometry[MESH_0], app->solver[MESH_0], app->config, MESH_0, 0,0,false);
      
      /* Compute eddy viscosity (needed by the flow solver) */
      
      app->solver[MESH_0][TURB_SOL]->Postprocessing(app->geometry[MESH_0], app->solver[MESH_0], app->config, MESH_0);      
    }
    
    
    /* Restrict flow and turb. solution to coarser meshes */
    
    for (unsigned short iMGLevel = 0; iMGLevel < app->config->GetnMGLevels(); iMGLevel++) {
      
      app->integration[FLOW_SOL]->SetRestricted_Solution(RUNTIME_FLOW_SYS, app->solver[iMGLevel][FLOW_SOL],
                             app->solver[iMGLevel+1][FLOW_SOL],
                             app->geometry[iMGLevel], app->geometry[iMGLevel+1], app->config);
        
      if (turbulent){
        app->integration[TURB_SOL]->SetRestricted_Solution(RUNTIME_TURB_SYS, app->solver[iMGLevel][TURB_SOL], app->solver[iMGLevel+1][TURB_SOL], app->geometry[iMGLevel], app->geometry[iMGLevel+1], app->config);
        
        app->integration[TURB_SOL]->SetRestricted_EddyVisc(RUNTIME_TURB_SYS, app->solver[iMGLevel][TURB_SOL], app->solver[iMGLevel+1][TURB_SOL], app->geometry[iMGLevel], app->geometry[iMGLevel+1], app->config);
        
      }
    }
    
    
    /* Set solution_time_n on all grids */
    
    for (unsigned short iMesh = 0; iMesh <= app->config->GetnMGLevels(); iMesh++) {
      app->integration[FLOW_SOL]->SetDualTime_Solver(app->geometry[iMesh], app->solver[iMesh][FLOW_SOL], app->config, iMesh);
    }
    
    if (turbulent){
      
      app->integration[TURB_SOL]->SetDualTime_Solver(app->geometry[MESH_0], app->solver[MESH_0][TURB_SOL], app->config, MESH_0);
      
    }

    /* Take the next time step to tstop */
    if (app->rank_x == MASTER_NODE) cout << app->rank_t << ": STEP to " << tstop << ", dt = " << deltat
                                         << " accur " << app->config->GetOrderMagResidual() << " " << app->config->GetMinLogResidual() << endl;
    app->driver->Run();
    
    
    if (app->config->GetBraid_CoarseGrid_Space()){ 
      if (level > 0){
        app->integration[FLOW_SOL]->SetProlongated_Solution(RUNTIME_FLOW_SYS, app->solver[app->config->GetFinestMesh()-1][FLOW_SOL],
            app->solver[app->config->GetFinestMesh()][FLOW_SOL],
            app->geometry[app->config->GetFinestMesh()-1], app->geometry[app->config->GetFinestMesh()],
            app->config);
        app->config->SetFinestMesh(FinestMesh);
        
      }
    }

    /* Check for SU2 convergence */
    if (!app->integration[FLOW_SOL]->GetConvergence() && (app->rank_x == MASTER_NODE)) {
        cout<<"ERROR: SU2 Solver didn't converge!? resid[0]: " << log10(SU2_TYPE::GetValue(app->solver[MESH_0][FLOW_SOL]->GetRes_RMS(0))) << endl;
        //exit(EXIT_FAILURE);
    }

    /* Update the Solution_n and solution_n1 for dual time stepping strategy */
    app->driver->Update();

    /* Grab the history values from SU2's master node. */
    if (app->rank_x == MASTER_NODE){
        u->Solution->Total_CD_n = SU2_TYPE::GetValue(app->solver[MESH_0][FLOW_SOL]->GetTotal_CD());
        u->Solution->Total_CL_n = SU2_TYPE::GetValue(app->solver[MESH_0][FLOW_SOL]->GetTotal_CL());
    }

    /* Grab the solution vectors from su2 for both time steps */
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
            if (app->BDF2) u->Solution->time_n1[iPoint][iVar] = SU2_TYPE::GetValue(app->solver[MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
      if (turbulent){
        for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
          u->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver[MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n()[iVar-nVar]);
          if (app->BDF2) u->Solution->time_n1[iPoint][iVar] = SU2_TYPE::GetValue(app->solver[MESH_0][TURB_SOL]->node[iPoint]->GetSolution_time_n1()[iVar-nVar]);
        }
        }
    }

    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nDim   = app->geometry[MESH_0]->GetnDim();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }
    

    /* Allocate memory for the primal braid vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar+nVar_Turb);

    /* Set the initial condition */
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar]  = app->initial_condition->time_n[iPoint][iVar];
            if(app->BDF2) u->Solution->time_n1[iPoint][iVar] = app->initial_condition->time_n1[iPoint][iVar];
        }
        if (turbulent){
          for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
              u->Solution->time_n[iPoint][iVar]  = app->initial_condition->time_n[iPoint][iVar];
              if(app->BDF2) u->Solution->time_n1[iPoint][iVar] = app->initial_condition->time_n1[iPoint][iVar];
          }
          
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}


int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nDim   = app->geometry[MESH_0]->GetnDim();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;

    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }
    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar+nVar_Turb);

    /* Copy the values from u to v */
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            if (app->BDF2) v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
        }
        if (turbulent){
          for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
              v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
              if (app->BDF2) v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
          }
        }
    }

    /* Set the pointer */
    *v_ptr = v;

    return 0;
}


int my_Free( braid_App app, braid_Vector u ){

    /* Delete the Solution vector */
    delete u->Solution;

    /* Delete the braid_Vector */
    delete u;

    return 0;
}



int my_Sum( braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }

    /* Compute the sum y = alpha x + beta y at time n and time n-1 */
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
                                               + beta  * y->Solution->time_n[iPoint][iVar];
            if (app->BDF2) y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
                                                              + beta  * y->Solution->time_n1[iPoint][iVar];
        }
        if (turbulent){
          for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
              y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
                                                 + beta  * y->Solution->time_n[iPoint][iVar];
              if (app->BDF2) y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
                                                                + beta  * y->Solution->time_n1[iPoint][iVar];
          }
        }
    }

    return 0;
}


int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }
    
    /* Compute l2norm of the solution list at time n and n1 */
    double norm = 0.0;
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->Solution->time_n[iPoint][iVar], 2);
            if (app->BDF2) norm += pow(u->Solution->time_n1[iPoint][iVar], 2);
        }
      if (turbulent){
        for (unsigned short iVar = nVar; iVar < nVar+nVar+nVar_Turb; iVar++){
          norm += pow(u->Solution->time_n[iPoint][iVar], 2);
          if (app->BDF2) norm += pow(u->Solution->time_n1[iPoint][iVar], 2);
        }
      }
    }

    /* Sum the norm from all spatial processors */
    double mynorm = norm;
    norm = 0.0;
    MPI_Allreduce(&mynorm, &norm, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    /* Set the pointer */
    *norm_ptr = sqrt(norm);

    return 0;
}

int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus ){

    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();

    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    su2double* cast = new su2double[nVar+nVar_Turb];
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }
    /* Get the current time-step number */
    double t;
    braid_AccessStatusGetT(astatus, &t);
    int iExtIter = (int) round( t / app->initialDT );

    if (iExtIter > 0){

        /* --- Write history and solution at time n ---*/

        app->config->SetExtIter(iExtIter);

        /* Trick SU2 with the current solution */
        for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
            for (unsigned short iVar = 0; iVar < nVar; iVar++){
                cast[iVar]  = u->Solution->time_n[iPoint][iVar];
            }
            app->solver[MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
            if (turbulent){
              for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
                  cast[iVar-nVar]  = u->Solution->time_n[iPoint][iVar];
              }
              app->solver[MESH_0][TURB_SOL]->node[iPoint]->SetSolution(cast);
            }
        }

        app->solver[MESH_0][FLOW_SOL]->Preprocessing(app->geometry[MESH_0], app->solver[MESH_0], app->config, MESH_0, 0,0,false);
        if (turbulent){
          app->solver[MESH_0][FLOW_SOL]->Postprocessing(app->geometry[MESH_0], app->solver[MESH_0], app->config, MESH_0);
          
        }

        /*--- Calculate the inviscid and viscous forces ---*/
        app->solver[MESH_0][FLOW_SOL]->Pressure_Forces(app->geometry[MESH_0], app->config);
        app->solver[MESH_0][FLOW_SOL]->Momentum_Forces(app->geometry[MESH_0], app->config);
        app->solver[MESH_0][FLOW_SOL]->Friction_Forces(app->geometry[MESH_0], app->config);


        /* Write to the history file */
        app->driver->Monitor(iExtIter);
        
        app->driver->Output(iExtIter);

        /* TODO: Also write time_n1 if BDF2!!  */
    }

    delete [] cast;

    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }
    /* Compute size of buffer */
    *size_ptr = nPoint * (nVar+nVar_Turb) * sizeof(double);
    if (app->BDF2) {
        *size_ptr = 2.0 * nPoint * (nVar+nVar_Turb) * sizeof(double);
    }


    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    unsigned long nPoint = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }

    /* Pack the buffer with current and previous time */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            /* Write Solution at current time to the buffer */
            dbuffer[ibuffer] = u->Solution->time_n[iPoint][iVar];
            ibuffer++;
            if (app->BDF2){
                /* Write Solution at previous time to the buffer */
                dbuffer[ibuffer] = u->Solution->time_n1[iPoint][iVar];
                ibuffer++;
            }
        }
        if (turbulent){
          for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
              /* Write Solution at current time to the buffer */
              dbuffer[ibuffer] = u->Solution->time_n[iPoint][iVar];
              ibuffer++;
              if (app->BDF2){
                  /* Write Solution at previous time to the buffer */
                  dbuffer[ibuffer] = u->Solution->time_n1[iPoint][iVar];
                  ibuffer++;
              }
            }
        }
    }

    /* Set the size of the buffer */
    int size = nPoint * (nVar+nVar_Turb) * sizeof(double);
    if (app->BDF2) size = 2.0 * size;
    braid_BufferStatusSetSize( bstatus, size );


    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    unsigned long nPoint  = app->geometry[MESH_0]->GetnPoint();
    unsigned short nVar   = app->solver[MESH_0][FLOW_SOL]->GetnVar();
    unsigned short nVar_Turb = 0;
    bool turbulent = app->config->GetKind_Turb_Model() != NONE;
    if (turbulent){
      nVar_Turb = app->solver[MESH_0][TURB_SOL]->GetnVar();
    }

    /* Get the buffer */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;

    /* Allocate memory for the new braid Vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    /* Unpack the buffer and write solution at current and previous time */
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        for (unsigned short iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar] = dbuffer[ibuffer];
            ibuffer++;
            if (app->BDF2){
                u->Solution->time_n1[iPoint][iVar] = dbuffer[ibuffer];
                ibuffer++;
            }
        }
        if (turbulent){
          for (unsigned short iVar = nVar; iVar < nVar+nVar_Turb; iVar++){
              u->Solution->time_n[iPoint][iVar] = dbuffer[ibuffer];
              ibuffer++;
              if (app->BDF2){
                  u->Solution->time_n1[iPoint][iVar] = dbuffer[ibuffer];
                  ibuffer++;
              }
            }
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}


/**
 * Conversion function which uses a format specifier for the string conversion.
 *
 * @param format        The format specifier like printf
 * @param list          The variable argument list for the format string
 *
 * @return  The output with the formated values
 */
std::string vformat(const char* format, va_list list) {
    const int bufferSize = 200;
    char buffer[bufferSize];

    // copy the list if we need to iterate through the variables again
    va_list listCpy;
    va_copy(listCpy, list);


    int outSize = vsnprintf(buffer, bufferSize, format, list);

    std::string result;
    if(outSize + 1 > bufferSize) {
        char* newBuffer = new char[outSize + 1];

        outSize = vsnprintf(newBuffer, outSize + 1, format, listCpy);

        result = newBuffer;

        delete [] newBuffer;
    } else {
        result = buffer;
    }

    // cleanup the copied list
    va_end (listCpy);

    return result;
}

/**
 * Conversion function which uses a format specifier for the string conversion.
 *
 * @param format        The format specifier like printf
 * @param ...           The values for the format string
 *
 * @return  The output with the formated values
 */
std::string format(const char* format, ...) {
    va_list list;
    va_start(list, format);
    std::string output = vformat(format, list);
    va_end(list);

    return output;
}

