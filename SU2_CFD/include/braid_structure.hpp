/*!
 * \file braid_structure.hpp
 * \brief Headers of structures and function for XBraid Integration
 *        The functions are in the <i>braid_wrapper.cpp</i> file.
 * \author S. Guenther
 *
 */

#pragma once

#include <braid.hpp>
#include "../../Common/include/mpi_structure.hpp"
#include "driver_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "transfer_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/interpolation_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include <memory>


/*!
 * \brief XBraid structure that holds additional information needed to carry out an unseady simulation step.
 */
typedef struct _braid_App_struct
{
  double tstart       = 0.0;  /* Begin of Time integration */
  double tstop        = 0.0;  /* End of Time integration */
  int    ntime        = 0;    /* Number of time steps */
  double initialDT    = 0.0;  /* Initial DeltaT */
  double initialstart = 0.0;  /* Initial starting time USED FOR TESTING ONLY */


  /* Information about communication */
  MPI_Comm comm;        /* global communicator */
  MPI_Comm comm_t;      /* temporal communicator */
  MPI_Comm comm_x;      /* spatial communicator */
  int su2rank;          /* Rank of current su2 processor */
  int braidrank;        /* Rank of current braid processor */


  /* Information about the braid grid */
  int ilower;              // smallest time index that belongs to this processor
  int iupper;              // largest time index that belongs to this processor


  /* Add the SU2 containers for SU2_CFD computations */
  /* TODO: Add only ZONE_0 ! (*_container[ZONE_0]) */
  CDriver *driver;
  CIteration **iteration_container;
  COutput *output;
  CIntegration ***integration_container;
  CGeometry ***geometry_container;
  CSolver ****solver_container;
  CNumerics *****numerics_container;
  CConfig **config_container;
  CSurfaceMovement **surface_movement;
  CVolumetricMovement **grid_movement;
  CFreeFormDefBox*** FFDBox;
  CInterpolator ***interpolator_container;
  CTransfer ***transfer_container;

  /* Output of history file */
  stringstream* history_stream;

  /* Information for optimization */
  double primal_norm    = 0.0;    // Norm of primal xBraid residual
  double Total_Cd_avg   = 0.0;    // Time-averaged objective funtion
  double Total_Cd_avg_b = 1.0;    // Seed for adjoint sensitivity computation
  double redgrad        = 0.0;    // Gradient of the costfunction wrt beta
  double redgrad_norm   = 0.0;    // Norm of the gradient
  int optimiter;                // Iteration number of outer optimization loop.


} my_App;

struct TwoStepSolution
{
    /* Dimensions of the solution lists */
    int nPoint;
    int nVar;
    /* Solution lists for all grid points and for all nVars*/
    double **time_n;    /*!<\brief List of solutions at time n for each point in space. */
    double **time_n1;   /*!<\brief List of solutions at time n-1 for each point in space. */

    /* Constructor */
    TwoStepSolution(int Point, int Var){
      nPoint = Point;
      nVar   = Var;
      /* Allocate memory for the solution lists */
      time_n  = new double*[nPoint];
      time_n1 = new double*[nPoint];
      for (int iPoint = 0; iPoint < nPoint; iPoint++){
          time_n[iPoint]  = new double[nVar];
          time_n1[iPoint] = new double[nVar];
      }
    }

    /* Destructor */
    ~TwoStepSolution(){
        for (int iPoint = 0; iPoint < nPoint; iPoint++){
          delete [] time_n[iPoint];
          delete [] time_n1[iPoint];
        }
        delete [] time_n;
        delete [] time_n1;
    }
};


/*!
 * \brief XBraid structure that defines a state vector at a certain time value and any information related to this vector which is needed to evolve the vector to the next time value, like mesh information.
  */
typedef struct _braid_Vector_struct
{

    /* Pointer to primal Solution lists for time steps n and n-1 */
    TwoStepSolution* Solution;

    /* Shared pointer to adjoint Solution lists at time n and n-1 */
    std::shared_ptr<TwoStepSolution> Solution_b;


    /* Flow solution coefficients for time n and time n1*/
    double Total_CLift_n, Total_CLift_n1;
    double Total_CDrag_n, Total_CDrag_n1;
    double Total_CSideForce_n, Total_CSideForce_n1;
    double Total_CEff_n, Total_CEff_n1;
    double Total_CMx_n, Total_CMx_n1;
    double Total_CMy_n, Total_CMy_n1;
    double Total_CMz_n, Total_CMz_n1;
    double Total_CFx_n, Total_CFx_n1;
    double Total_CFy_n, Total_CFy_n1;
    double Total_CFz_n, Total_CFz_n1;

    /* Flow residual */
    double residual_dens_n;
    double residual_dens_n1;

    /* Constructor for primal Solution list: Called by my_init */
    /* Constructor for adjoint Solution list: Called in my_init */
    /* Destructor for primal Solution list: Called by my_Free */
    /* Destructor for adjoint Solution list: Called automatically by shared pointer, if no pointer to the struct is left */

} my_Vector;

/*!
 * \brief Creates the braidTape
 */
void setupTapeData();

/*!
 * \brief Create a copy v of a given vector u
 * \param braid_App app
 * \param braid_Vector u that is to be copied
 * \return braid_Vector copy
 */
braid_Vector deep_copy( braid_App app, braid_Vector u );


/*!
 * \brief This function tells XBraid how to take a time step. It advances the vector u from tstart to tstop.
*/
int my_Step( braid_App        app,
             braid_Vector     ustop,
             braid_Vector     fstop,
             braid_Vector     u,
             braid_StepStatus status );

/*!
 *\brief Tells XBraid, how to initialize a vector at time t
 */
int my_Init( braid_App app, double t, braid_Vector *u_ptr );

/*!
 *\brief Tells XBraid, how to clone a vector u into a new vector, *v_ptr
 */
int my_Clone( braid_App app, braid_Vector u, braid_Vector *v_ptr );

/*!
 *\brief Tells XBraid, how to free a vector
 */
int my_Free( braid_App app, braid_Vector u );

/*!
 *\brief Tells XBraid, how to sum two vectors (y = ax + by)
 */
int my_Sum( braid_App app, double alpha, braid_Vector x, double beta,
    braid_Vector y );

/*!
 *\brief Tells XBraid, how to take the norm of a braid_Vector
 */
int my_SpatialNorm( braid_App app, braid_Vector u, double *norm_ptr );

/*!
 *\brief Allows the user access to XBraid and the current solution vector at time t.
 */
int my_Access( braid_App app, braid_Vector u, braid_AccessStatus astatus );

/*!
 *\brief XBraid function that computes the upper bound of the size of a solution vector.
 */
int my_BufSize ( braid_App app, int *size_ptr, braid_BufferStatus bstatus );

/*!
 *\brief XBraid function that packs a vector into a void * buffer for MPI communication
 */
int my_BufPack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus);

/*!
 *\brief XBraid function that unpacks a void * buffer into a vector
 */
int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus  );

enum class BraidCall_t {
  PHI       = 1,
  INIT      = 2,
  CLONE     = 3,
  FREE      = 4,
  SUM       = 5,
  BUFPACK   = 6,
  BUFUNPACK = 7,
  ACCESS    = 8,
};


struct BraidAction_t {
    BraidCall_t       braidCall;        /* The type of xBraid call */
    std::vector<int>  inIndex;          /* The indicees input xBraid vectors */
    int               outIndex;         /* The index of output xBraid vector */
    double            deltat;           /* The time step size that was used in the primal Step function*/
    braid_StepStatus  StepStatus;        /* The status of xBraid for Phi Calls */
    double            sum_alpha;        /* First coefficient of the sum */
    double            sum_beta;         /* Second coefficient of the sum */
    double            time;         /* Current number of time step, needed for access_adjoint */
    int               send_recv_rank;   /* Processor rank of the sender / receiver */
    int               optimiter;       /* Iteration number of xBraid */
    int               myid;             /* Processors id */

    /* Constructor */
    BraidAction_t() : inIndex() {
    }
};

struct BraidTape_t {
    std::vector<my_Vector*>     primal;   /* Intermediate primal braid vectors that are used in the nonlinear operations */
    std::vector<BraidAction_t>  action;   /* Actions during one xBraid iteration */

    std::vector<std::shared_ptr<TwoStepSolution>> adjoint; /* Intermediate adjoint braid vectors */

    /* Two more tapes that store pointers to xbraid input and output variables in each iteration */
   std::vector<std::shared_ptr<TwoStepSolution>> braid_input;
   std::vector<std::shared_ptr<TwoStepSolution>> braid_output;

    /* Constructor */
    BraidTape_t() : primal(), action(), adjoint() {}
//     in_Adjoint(),
//     out_Adjoint() {
//    }
};

// Define the tapes as extern
extern BraidTape_t* braidTape;

/*!
 * \brief Creates the braidTape
 */
void setupTapeData();

/*!
 * \brief Evaluates the Braid ActionTape in reverse order and calls adjoint actions
 */
void evalAdjointAction( braid_App app, BraidTape_t* braidTape);


/* Adjoint Action Calls */
void my_Step_adjoint( BraidAction_t &action, braid_App app );
void my_Access_adjoint( BraidAction_t &action , braid_App app );
void my_Sum_adjoint( BraidAction_t &action, braid_App app );
void my_Clone_adjoint( BraidAction_t &action, braid_App app );
void my_BufPack_adjoint( BraidAction_t &action, braid_App app );
void my_BufUnPack_adjoint( BraidAction_t &action, braid_App app );
