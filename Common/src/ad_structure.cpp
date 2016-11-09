/*!
 * \file dataype_structure.cpp
 * \brief Main subroutines for the datatype structures.
 * \author T. Albring
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/datatype_structure.hpp"

namespace AD {
#ifdef CODI_REVERSE_TYPE
  /*--- Initialization of the global variables ---*/

  int adjointVectorPosition = 0;

  std::vector<unsigned int> inputValues;
  std::vector<unsigned int> localInputValues;
  std::vector<su2double*> localOutputValues;

  codi::ChunkTape<double, int>& globalTape = codi::RealReverse::getGlobalTape();
  codi::ChunkTape<double, int>::Position StartPosition, EndPosition;

  bool Status = false;
  bool PreaccActive = false;

  void EndPreacc(){

    if(PreaccActive){
      unsigned short iVarOut, iVarIn;
      unsigned short nVarOut, nVarIn;
      int index_out, index_in;

      // nVarOut = localOutputValues.size();
      // nVarIn  = localInputValues.size();
      //
      // /*--- Store the current position of the tape ---*/
      //
      // EndPosition = globalTape.getPosition();
      //
      // /*--- Allocate local memory on the stack (does not need to be deleted at the end of the routine!) ---*/
      //
      // double* local_jacobi     = (double*)alloca(sizeof(double)*(nVarOut*nVarIn));
      // unsigned short* nNonzero = (unsigned short*)alloca(sizeof(unsigned short)*nVarOut);
      //
      // /*--- Compute the local Jacobi matrix of the code between the start and end position
      //  * using the inputs and outputs declared with StartPreacc(...)/EndPreacc(...) ---*/
      //
      // for (iVarOut = 0; iVarOut < nVarOut; iVarOut++){
      //   nNonzero[iVarOut] = 0;
      //   index_out = localOutputValues[iVarOut]->getGradientData();
      //
      //   globalTape.setGradient(index_out, 1.0);
      //   globalTape.evaluate(EndPosition, StartPosition);
      //
      //   for (iVarIn= 0; iVarIn < nVarIn; iVarIn++){
      //     index_in =  localInputValues[iVarIn];
      //     local_jacobi[iVarOut*nVarIn+iVarIn] = globalTape.getGradient(index_in);
      //     if (local_jacobi[iVarOut*nVarIn+iVarIn] != 0.0){
      //       nNonzero[iVarOut]++;
      //     }
      //     globalTape.setGradient(index_in, 0.0);
      //   }
      //   globalTape.setGradient(index_out, 0.0);
      //   globalTape.clearAdjoints(EndPosition, StartPosition);
      // }
      //
      // /*--- Reset the tape to the starting position (to reuse the part of the tape) ---*/
      //
      // if (nVarOut > 0){
      //   globalTape.reset(StartPosition);
      // }
      //
      // /*--- For each output create a statement on the tape and push the corresponding Jacobi entries.
      //  * Note that the output variables need a new index since we did a reset of the tape section. ---*/
      //
      // for (iVarOut = 0; iVarOut < nVarOut; iVarOut++){
      //   index_out = 0;
      //   if (nNonzero[iVarOut] != 0){
      //     globalTape.store(index_out, nNonzero[iVarOut]);
      //     for (iVarIn = 0; iVarIn < nVarIn; iVarIn++){
      //       index_in =  localInputValues[iVarIn];
      //      globalTape.pushJacobi(local_jacobi[iVarOut*nVarIn+iVarIn],
      //          local_jacobi[iVarOut*nVarIn+iVarIn], local_jacobi[iVarOut*nVarIn+iVarIn], index_in);
      //     }
      //   }
      //   localOutputValues[iVarOut]->getGradientData() = index_out;
      // }
      //
      // /* --- Clear local vectors and reset indicator ---*/
      //
      // localInputValues.clear();
      // localOutputValues.clear();
      //
      PreaccActive = false;
    }
  }
#endif
}
