/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_CFD.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  
  unsigned short nZone, nDim;
  char config_file_name[MAX_STRING_SIZE];
  bool fsi, turbo, xbraid;
  
  /*--- MPI initialization, and buffer setting ---*/

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
#ifdef HAVE_MPI
  int  buffsize;
  char *buffptr;
  SU2_MPI::Init(&argc, &argv);
  MPI_Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  SU2_MPI::comm = MPI_COMM_WORLD;
  SU2_MPI::comm_x = MPI_COMM_WORLD;
  SU2_MPI::comm_t = MPI_COMM_WORLD;
#else
  SU2_Comm MPICommunicator(0);
#endif
  
  /*--- Create a pointer to the main SU2 Driver ---*/
  
  CDriver *driver = NULL;

  /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
   file is specified, default.cfg is used) ---*/

  if (argc == 2) { strcpy(config_file_name, argv[1]); }
  else { strcpy(config_file_name, "default.cfg"); }

  /*--- Read the name and format of the input mesh file to get from the mesh
   file the number of zones and dimensions from the numerical grid (required
   for variables allocation)  ---*/

  CConfig *config = NULL;
  config = new CConfig(config_file_name, SU2_CFD);

  nZone  = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
  nDim   = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
  fsi    = config->GetFSI_Simulation();
  turbo  = config->GetBoolTurbomachinery();
  xbraid = config->GetBraid_Run();



    /* --- Preprocess the processor grid --- */

    if ( config->GetBraid_Run() ){
      if ( size % config->GetBraid_NProc_Time() != 0 ){
        cout << "\n\nError: px*pt does not equal the number of processors!\n\n";
        exit(EXIT_FAILURE);
      } else {
        /* Split communicators for the time and space dimensions */
        int px = size / config->GetBraid_NProc_Time();
        braid_SplitCommworld(&(SU2_MPI::comm), px, &(SU2_MPI::comm_x), &(SU2_MPI::comm_t));
        /* Pass the spatial communicator to SU2 */
//        SU2_MPI::comm_x = comm_x;
        /* Get the rank and size of braid and su2 processors */
//        MPI_Comm_size(comm_t, &braidsize);
//        MPI_Comm_size(comm_x, &su2size);
//        MPI_Comm_rank(comm_t, &braidrank);
//        MPI_Comm_rank(comm_x, &su2rank);
    }
}




  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem
   and perform all the preprocessing. ---*/

  if ( (config->GetKind_Solver() == FEM_ELASTICITY ||
        config->GetKind_Solver() == POISSON_EQUATION ||
        config->GetKind_Solver() == WAVE_EQUATION ||
        config->GetKind_Solver() == HEAT_EQUATION) ) {

    /*--- Single zone problem: instantiate the single zone driver class. ---*/
    
    if (nZone > 1 ) {
      cout << "The required solver doesn't support multizone simulations" << endl; 
      exit(EXIT_FAILURE);
    }
    
    driver = new CGeneralDriver(config_file_name, nZone, nDim, MPICommunicator);

  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {

    /*--- Harmonic balance problem: instantiate the Harmonic Balance driver class. ---*/

    driver = new CHBDriver(config_file_name, nZone, nDim, MPICommunicator);

  } else if ((nZone == 2) && fsi) {

    /*--- FSI problem: instantiate the FSI driver class. ---*/

    driver = new CFSIDriver(config_file_name, nZone, nDim, MPICommunicator);

  } else {

    /*--- Multi-zone problem: instantiate the multi-zone driver class by default
    or a specialized driver class for a particular multi-physics problem. ---*/

    if (config->GetDiscrete_Adjoint()) {

      if (turbo) {

        driver = new CDiscAdjTurbomachineryDriver(config_file_name, nZone, nDim, MPICommunicator);

      } else {

        driver = new CDiscAdjFluidDriver(config_file_name, nZone, nDim, MPICommunicator);
        
      }

    } else if (turbo) {

      driver = new CTurbomachineryDriver(config_file_name, nZone, nDim, MPICommunicator);

    } else {

      /*--- Instantiate the class for external aerodynamics ---*/

      driver = new CFluidDriver(config_file_name, nZone, nDim, MPICommunicator);
      
    }
    
  }

  delete config;
  config = NULL;


  /*--- Launch the main external loop of the solver ---*/
  
  if( !xbraid ) {

      /* Time-serial run */

      driver->StartSolver();

  } else {


      /* Time-parallel XBraid iteration */

      driver->StartXBraidSolver();

//      for (int optimiter = 0; optimiter < config_container[ZONE_0]->GetBraid_Max_Iter(); optimiter++){

//        /* Reset the app */
//        app->Total_Cd_avg   = 0.0;
//        app->Total_Cd_avg_b = 1.0;
//        app->optimiter      = optimiter;


//        /* Clear the action tape */
//        braidTape->action.clear();

//        /* --- Primal xBraid computation ---*/

//        /* Run one primal xBraid iteration */
//        braid_Drive(driver->xbraidcore);


//        /* Get the primal xBraid residuum */
//        _braid_GetRNorm(xbraidcore, -1, &app->primal_norm);

//        /* Compute the time-average of CDrag */
//        double MyTotalAvg = app->Total_Cd_avg;
//        app->Total_Cd_avg = 0.0;
//        MPI_Allreduce(&MyTotalAvg, &app->Total_Cd_avg, 1, MPI_DOUBLE, MPI_SUM, comm_t);
//        app->Total_Cd_avg = 1.0/(app->ntime * 2) * app->Total_Cd_avg;
//        cout<< format("Total_Cd_avg %1.14e\n", app->Total_Cd_avg);


//        /* --- Adjoint sensitivity computation --- */
//        cout<< "ADJONT\n";

//        /* Reset the reduced gradient */
//        for (int iPoint = 0; iPoint < nPoint; iPoint++){
//          for (int iDim = 0; iDim < nDim; iDim++){
//            app->redgrad[iPoint][iDim] = 0.0;
//          }
//        }

//        /* Adjoint of computing the time-average. */
//         app->Total_Cd_avg_b = 1.0/(app->ntime * 2 ) * app->Total_Cd_avg_b;

//        /* Evaluate the Action tape in reverse order. */
//        evalAdjointAction(app, braidTape);

//        /* Compute adjoint residuum */
//        double my_norm = 0.0;
//        for (int i = 0; i < app->ncpoints; i++)
//        {
//          /* TODO: COMPUTE THE NORM! */
//          // if (myid != 0 || i !=0 )
//          // {
//            // my_norm += pow(getValue(braidTape->braid_input_b[i]->y) - getValue(app->optim->adjoint[i]->y), 2);
//          // }
//        }
//        MPI_Allreduce(&my_norm, &app->adjoint_norm, 1, MPI_DOUBLE, MPI_SUM, comm);
//        app->adjoint_norm = sqrt(app->adjoint_norm);
//        // if (optimiter == 0) app->optim_adjoint_norm0 = app_optim->adjoint_norm;

//        /* Compute the reduced gradient */
//        // double MyRedGrad = app->redgrad;
//        // app->redgrad= 0.0;
//        // MPI_Allreduce(&MyRedGrad, &app->redgrad, 1, MPI_DOUBLE, MPI_SUM, comm);

//        /* Store the adjoints into the Optim structure */
//        for (int i=0; i < app->ncpoints; i++)
//        {
//          /* TODO: Store the adjoint in optimadjoint
//          // app->optim->adjoint[i]->y = braidTape->in_Adjoint[i]->y;
//          /* Delete the pointer */
//          braidTape->braid_input_b[i].reset();
//        }

//        /* Move pointers from braid_output_b to braid_input_b */
//        // Because braid output of current iteration is braid input of next iteration
//        for (int i=0; i<app->ncpoints; i++)
//        {
//          braidTape->braid_input_b[i] = braidTape->braid_output_b[i];
//          braidTape->braid_output_b[i].reset();
//        }

//        /* Output */
//        if (rank == MASTER_NODE){
//          cout<<format(" || r_%d || = %1.14e  CD_avg = %1.14e\n", optimiter, app->primal_norm, app->Total_Cd_avg);
//        }

//        /* Stopping criterion */
//        if (app->primal_norm < app->config_container[ZONE_0]->GetBraid_Tol()){
//            cout<< format("\n XBraid has converged! primal res = %1.14e \n\n", app->primal_norm);
//            break;
//        }

//      } // END OF OPTIMIZATION LOOP

//      /* Print some statistics */
//      braid_PrintStats(xbraidcore);

//    } else {
//      std::cout<<format("\n\nTurn warm_restart option on for One-Shot!!\n\n");
//      return -1;
//    }


//  /* Print the sensitivity of a surface point */
//  /* For finite differencing only!! */
//  for (int iMarker = 0; iMarker < app->geometry_container[ZONE_0][MESH_0]->GetnMarker(); iMarker++){
//    if(app->config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == EULER_WALL
//        || app->config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
//        || app->config_container[ZONE_0]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL){
//      int iPoint_vertex0 = app->geometry_container[ZONE_0][MESH_0]->vertex[iMarker][0]->GetNode();
//      cout<< format("grad surface %d %1.14e\n", iPoint_vertex0, app->redgrad[iPoint_vertex0][0]);
//    }
//  }


////    std::ofstream out;
////    ParallelFileIO::startFileWrite(out, "history_test.dat", braidrank, braidsize, 42, comm_t);
////    cout << (*app->history_stream).str();
////    out << (*app->history_stream).str();
////    ParallelFileIO::endFileWrite(out, braidrank, braidsize, 42, comm_t);


//    // Finalize XBraid
//    braid_Destroy(xbraidcore);

  }



  /*--- Postprocess all the containers, close history file, exit SU2 ---*/
  
  driver->Postprocessing();

  if (driver != NULL) delete driver;
  driver = NULL;

  /*--- Finalize MPI parallelization ---*/

#ifdef HAVE_MPI
  MPI_Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
