/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author F. Palacios, T. Economon
 * \version 6.0.1 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
  bool fsi, turbo, zone_specific, periodic = false, xbraid;
  
  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
  int  buffsize;
  char *buffptr;
  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
  SU2_MPI::Comm comm_su2, comm_braid;
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

  nZone    = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
  nDim     = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
  fsi      = config->GetFSI_Simulation();
  turbo    = config->GetBoolTurbomachinery();
  periodic = CConfig::GetPeriodic(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
  zone_specific = config->GetBoolZoneSpecific();
  xbraid   = config->GetBraid_Run();



  /* --- Preprocess the processor grid --- */

  if ( config->GetBraid_Run() ){
    if ( SU2_MPI::GetGlobalSize() % config->GetBraid_NProc_Time() != 0 ){
      SU2_MPI::Error(" XBraid: npx * npt does not equal to the number of processes used! ", CURRENT_FUNCTION);
    } else {
      /* Split communicators for the time and space dimensions */
      int px = SU2_MPI::GetGlobalSize()/ config->GetBraid_NProc_Time();
      SU2_MPI::Comm comm = MPI_COMM_WORLD;
      braid_SplitCommworld(&(comm), px, &comm_su2, &comm_braid);
      SU2_MPI::SetComm(comm_su2);
      SU2_MPI::SetTimeComm(comm_braid);
      MPICommunicator = comm_su2;
    }
  }


  //csg: TESTING
  int rank_su2, size_su2, rank_braid, size_braid, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rank_su2   = SU2_MPI::GetRank();
  rank_braid = SU2_MPI::GetGlobalRank();
  cout << rank << ": " << rank_braid << " "<< rank_su2 << endl;


  /*--- First, given the basic information about the number of zones and the
   solver types from the config, instantiate the appropriate driver for the problem
   and perform all the preprocessing. ---*/

  if ( (config->GetKind_Solver() == FEM_ELASTICITY ||
        config->GetKind_Solver() == DISC_ADJ_FEM ||
        config->GetKind_Solver() == POISSON_EQUATION ||
        config->GetKind_Solver() == WAVE_EQUATION ||
        config->GetKind_Solver() == HEAT_EQUATION) ) {

    /*--- Single zone problem: instantiate the single zone driver class. ---*/
    
    if (nZone > 1 ) {
      SU2_MPI::Error("The required solver doesn't support multizone simulations", CURRENT_FUNCTION);
    }
    
    driver = new CGeneralDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);

  } else if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {

    /*--- Harmonic balance problem: instantiate the Harmonic Balance driver class. ---*/

    driver = new CHBDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);

  } else if ((nZone == 2) && fsi) {

    bool stat_fsi = ((config->GetDynamic_Analysis() == STATIC) && (config->GetUnsteady_Simulation() == STEADY));
    bool disc_adj_fsi = (config->GetDiscrete_Adjoint());

    /*--- If the problem is a discrete adjoint FSI problem ---*/
    if (disc_adj_fsi) {
      if (stat_fsi) {
        driver = new CDiscAdjFSIDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);
      }
      else {
        SU2_MPI::Error("WARNING: There is no discrete adjoint implementation for dynamic FSI. ", CURRENT_FUNCTION);
      }
    }
    /*--- If the problem is a direct FSI problem ---*/
    else{
      driver = new CFSIDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);
    }

  } else if (zone_specific) {
    driver = new CMultiphysicsZonalDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);
  } else {

    /*--- Multi-zone problem: instantiate the multi-zone driver class by default
    or a specialized driver class for a particular multi-physics problem. ---*/

    if (config->GetDiscrete_Adjoint()) {

      if (turbo) {

        driver = new CDiscAdjTurbomachineryDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);

      } else {

        driver = new CDiscAdjFluidDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);
        
      }

    } else if (turbo) {

      driver = new CTurbomachineryDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);

    } else {

      /*--- Instantiate the class for external aerodynamics ---*/

      driver = new CFluidDriver(config_file_name, nZone, nDim, periodic, MPICommunicator);
      
    }
    
  }

  delete config;
  config = NULL;


  /*--- Launch the main external loop of the solver ---*/
  
  // driver->StartSolver();

  /*--- Postprocess all the containers, close history file, exit SU2 ---*/
  
  driver->Postprocessing();

  if (driver != NULL) delete driver;
  driver = NULL;

  /*--- Finalize MPI parallelization ---*/

#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif
  
  return EXIT_SUCCESS;
  
}
