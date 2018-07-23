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
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        /* Copy the values from u to v */
        /* TODO: Set the values with the SetSolution(...), Set_solution_time_n(...) and Set_Solution_timen1(...) */
        for (int iVar = 0; iVar < nVar; iVar++){
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

    double  tstart, tstop, deltat;
    int     nPoint, nVar;
    int     iExtIter;

    /* Grab variables from the app */
    nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  
    /* Set the time-step size and time-step index */
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
    if (app->BDF2) {
        deltat = (tstop - tstart) / 2.0;
    } else {
        deltat = (tstop - tstart);
    }
    iExtIter = (int) round( tstop  / app->initialDT );
    app->config_container[ZONE_0]->SetDelta_UnstTimeND( deltat );
    app->config_container[ZONE_0]->SetExtIter(iExtIter); //TODO: Check if BDF2

    /* Set the solution vectors (Solution, Solution_time_n and Solution_time_n1*/
    su2double *cast_n, *cast_n1;
    cast_n = new su2double[nVar];
    if(app->BDF2) cast_n1 = new su2double[nVar];
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            cast_n[iVar]  = u->Solution->time_n[iPoint][iVar];
            if(app->BDF2) cast_n1[iVar] = u->Solution->time_n1[iPoint][iVar];
        }
        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast_n);
        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(cast_n);
        if (app->BDF2) app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(cast_n1);
    }
    delete [] cast_n;
    if(app->BDF2) delete [] cast_n1;

    /* Interpolate the solution_rime_n and solution_time_n1 to all meshes */
    su2double *Solution_n, *Solution_Fine_n, *Solution_n1, *Solution_Fine_n1, Area_Parent, Area_Children;
    unsigned Point_Fine;
    Solution_n = new su2double[nVar];
    Solution_n1 = new su2double[nVar];
    for (unsigned short iMesh = 1; iMesh <= app->config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
        for (int iPoint = 0; iPoint < app->geometry_container[ZONE_0][INST_0][iMesh]->GetnPoint(); iPoint++) {
            Area_Parent = app->geometry_container[ZONE_0][INST_0][iMesh]->node[iPoint]->GetVolume();
            for (int iVar = 0; iVar < nVar; iVar++) {
                Solution_n[iVar]  = 0.0;
                Solution_n1[iVar] = 0.0;
            }
            for (int iChildren = 0; iChildren < app->geometry_container[ZONE_0][INST_0][iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
                Point_Fine = app->geometry_container[ZONE_0][INST_0][iMesh]->node[iPoint]->GetChildren_CV(iChildren);
                Area_Children = app->geometry_container[ZONE_0][INST_0][iMesh-1]->node[Point_Fine]->GetVolume();
                Solution_Fine_n  = app->solver_container[ZONE_0][INST_0][iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution_time_n();
                Solution_Fine_n1 = app->solver_container[ZONE_0][INST_0][iMesh-1][FLOW_SOL]->node[Point_Fine]->GetSolution_time_n1();
                for (int iVar = 0; iVar < nVar; iVar++) {
                    Solution_n[iVar]  += Solution_Fine_n[iVar]*Area_Children/Area_Parent;
                    Solution_n1[iVar] += Solution_Fine_n1[iVar]*Area_Children/Area_Parent;
                }
            }
            app->solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->node[iPoint]->SetSolution_time_n(Solution_n);
            app->solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->node[iPoint]->SetSolution_time_n1(Solution_n1);
        }
        app->solver_container[ZONE_0][INST_0][iMesh][FLOW_SOL]->Set_MPI_Solution(app->geometry_container[ZONE_0][INST_0][iMesh], app->config_container[ZONE_0]);
    }
    delete [] Solution_n;
    delete [] Solution_n1;


    if (app->BDF2){

        /* Take the first time step to tstart + deltat */
        // if (app->rank_x == MASTER_NODE) cout << app->rank_t << ": STEP from " << tstart << " to " << tstart + deltat <<  endl;
        app->driver->Run();

        /* Check for SU2 convergence */
        if (!app->integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence()) {
            cout<<format("ERROR: SU2 Solver didn't converge!? resid[0]: %1.14e\n", SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)));
            exit(EXIT_FAILURE);
        }

        /* Update the Solution_n and solution_n1 for dual time stepping strategy */
        app->driver->Update();

        /* Grab drag and lift coefficient from SU2's master node. */
        if (app->rank_x == MASTER_NODE){
            u->Solution->Total_CD_n1 = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD());
            u->Solution->Total_CL_n1 = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL());
        }

        /* Set the next iExtIter */
        iExtIter++;
        app->config_container[ZONE_0]->SetExtIter(iExtIter);
    }


    /* Take the next time step to tstop */
    if (app->rank_x == MASTER_NODE) cout << app->rank_t << ": STEP to " << tstop << ", dt = " << deltat <<  endl;
    app->driver->Run();

    /* Check for SU2 convergence */
    if (!app->integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence()) {
        cout<<format("ERROR: SU2 Solver didn't converge!? resid[0]: %1.14e\n", SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0)));
        //exit(EXIT_FAILURE);
    }

    /* Update the Solution_n and solution_n1 for dual time stepping strategy */
    app->driver->Update();

    /* Grab the history values from SU2's master node. */
    if (app->rank_x == MASTER_NODE){
        u->Solution->Total_CD_n = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD());
        u->Solution->Total_CL_n = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL());
    }

    /* Grab the solution vectors from su2 for both time steps */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar]  = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n()[iVar]);
            if (app->BDF2) u->Solution->time_n1[iPoint][iVar] = SU2_TYPE::GetValue(app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1()[iVar]);
        }
    }

    return 0;
}


int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nDim   = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Allocate memory for the primal braid vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    /* Set the initial condition */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar]  = app->initial_condition->time_n[iPoint][iVar];
            if(app->BDF2) u->Solution->time_n1[iPoint][iVar] = app->initial_condition->time_n1[iPoint][iVar];
        }
    }

    /* Set the pointer */
    *u_ptr = u;

    return 0;
}


int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nDim   = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnDim();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();


    /* Allocate memory for the new copy v */
    my_Vector* v = new my_Vector;
    v->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    /* Copy the values from u to v */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            v->Solution->time_n[iPoint][iVar]  = u->Solution->time_n[iPoint][iVar];
            if (app->BDF2) v->Solution->time_n1[iPoint][iVar] = u->Solution->time_n1[iPoint][iVar];
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
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();


    /* Compute the sum y = alpha x + beta y at time n and time n-1 */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            y->Solution->time_n[iPoint][iVar]  = alpha * x->Solution->time_n[iPoint][iVar]
                                               + beta  * y->Solution->time_n[iPoint][iVar];
            if (app->BDF2) y->Solution->time_n1[iPoint][iVar] = alpha * x->Solution->time_n1[iPoint][iVar]
                                                              + beta  * y->Solution->time_n1[iPoint][iVar];
        }
    }

    return 0;
}


int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute l2norm of the solution list at time n and n1 */
    double norm = 0.0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            norm += pow(u->Solution->time_n[iPoint][iVar], 2);
            if (app->BDF2) norm += pow(u->Solution->time_n1[iPoint][iVar], 2);
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

    int nPoint      = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar        = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
    su2double* cast = new su2double[nVar];

    /* Get the current time-step number */
    double t;
    braid_AccessStatusGetT(astatus, &t);
    int iExtIter = (int) round( t / app->initialDT );

    if (iExtIter > 0){

        /* --- Write history and solution at time n ---*/

        app->config_container[ZONE_0]->SetExtIter(iExtIter);

        /* Trick SU2 with the current solution */
        for (int iPoint = 0; iPoint < nPoint; iPoint++){
            for (int iVar = 0; iVar < nVar; iVar++){
                cast[iVar]  = u->Solution->time_n[iPoint][iVar];
            }
            app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(cast);
        }

        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Preprocessing(app->geometry_container[ZONE_0][INST_0][MESH_0], app->solver_container[ZONE_0][INST_0][MESH_0], app->config_container[ZONE_0], MESH_0, 0,0,false);

//        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->Set_MPI_Solution(app->geometry_container[ZONE_0][MESH_0],app->config_container[ZONE_0]);



        /* Compute the primitive Variables from the conservative ones */
//        app->solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetPrimitive_Variables(app->solver_container[ZONE_0][MESH_0], app->config_container[ZONE_0], false);

        /*--- Calculate the inviscid and viscous forces ---*/
        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Pressure_Forces(app->geometry_container[ZONE_0][INST_0][MESH_0], app->config_container[ZONE_0]);
        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Momentum_Forces(app->geometry_container[ZONE_0][INST_0][MESH_0], app->config_container[ZONE_0]);
        app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->Friction_Forces(app->geometry_container[ZONE_0][INST_0][MESH_0], app->config_container[ZONE_0]);


        /* Write to the history file */
        app->driver->Monitor(iExtIter);
  
        unsigned short nInst[] = {1};

        /* Write the solution files */
        app->output->SetResult_Files_Parallel(app->solver_container, app->geometry_container,
                                              app->config_container, iExtIter, 1, nInst);

        /* TODO: Also write time_n1 if BDF2!!  */
    }

    delete [] cast;

    return 0;
}

int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Compute size of buffer */
    *size_ptr = nPoint * nVar * sizeof(double);
    if (app->BDF2) {
        *size_ptr = 2.0 * nPoint * nVar * sizeof(double);
    }


    return 0;
}

int my_BufPack( braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar   = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();


    /* Pack the buffer with current and previous time */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
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

    /* Set the size of the buffer */
    int size = nPoint * nVar * sizeof(double);
    if (app->BDF2) size = 2.0 * size;
    braid_BufferStatusSetSize( bstatus, size );


    return 0;
}

int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  ){

    /* Grab variables from the app */
    int nPoint  = app->geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();
    int nVar    = app->solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();

    /* Get the buffer */
    double *dbuffer = (double*)buffer;
    int ibuffer = 0;

    /* Allocate memory for the new braid Vector */
    my_Vector* u = new my_Vector;
    u->Solution  = new TwoStepSolution(app->BDF2, nPoint, nVar);

    /* Unpack the buffer and write solution at current and previous time */
    for (int iPoint = 0; iPoint < nPoint; iPoint++){
        for (int iVar = 0; iVar < nVar; iVar++){
            u->Solution->time_n[iPoint][iVar] = dbuffer[ibuffer];
            ibuffer++;
            if (app->BDF2){
                u->Solution->time_n1[iPoint][iVar] = dbuffer[ibuffer];
                ibuffer++;
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

