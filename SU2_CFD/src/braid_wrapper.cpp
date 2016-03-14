/*!
 * \file braid_wrapper.cpp
 * \brief Functions for XBraid integration.
 *
 * \author S. Guenther
 *
 */

 #include <braid.hpp>

int my_Phi( braid_App app, braid_Vector u, braid_PhiStatus status ){

  return 0;
}

int my_Init( braid_App app, double t, braid_Vector *u_ptr ){

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
