/*!
 * \file braid_structure.hpp
 * \brief Headers of structures and function for XBraid Integration
 *        The functions are in the <i>braid_wrapper.cpp</i> file.
 * \author S. Guenther
 *
 */

#pragma once

#include <braid.hpp>

/*!
 * \brief XBraid structure that holds additional information needed to carry out an unseady simulation step.
 */
typedef struct _braid_App_struct
{
  su2double tstart;     /* Begin of Time integration */
  su2double tstop;      /* End of Time integration */
  int ntime;            /* Number of time steps */

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

} my_App;

/*!
 * \brief XBraid structure that defines a state vector at a certain time value and any information related to this vector which is needed to evolve the vector to the next time value, like mesh information.
  */
typedef struct _braid_Vector_struct
{
  CVariable** node;	    /*!< \brief Vector which defines the flow variables for each problem. */

} my_Vector;

/*!
 * \brief This function tells XBraid how to take a time step. It advances the vector u from tstart to tstop.
*/
int my_Phi( braid_App app, braid_Vector u, braid_PhiStatus status );

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
int my_BufSize ( braid_App app, int *size_ptr );

/*!
 *\brief XBraid function that packs a vector into a void * buffer for MPI communication
 */
int my_BufPack( braid_App app, braid_Vector u, void *buffer,
                braid_Int *size_ptr );

/*!
 *\brief XBraid function that unpacks a void * buffer into a vector
 */
int my_BufUnpack( braid_App app, void *buffer, braid_Vector *u_ptr );
