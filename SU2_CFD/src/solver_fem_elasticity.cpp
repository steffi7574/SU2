/*!
 * \file solution_template.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
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

#include "../include/solver_structure.hpp"

CFEM_ElasticitySolver::CFEM_ElasticitySolver(void) : CSolver() {

	nElement = 0;
	nDim = 0;
	nMarker = 0;

	nPoint = 0;
	nPointDomain = 0;

	Total_CFEA = 0.0;
	WAitken_Dyn = 0.0;
	WAitken_Dyn_tn1 = 0.0;

	element_container = NULL;
	node = NULL;

	GradN_X = NULL;
	GradN_x = NULL;

	Jacobian_c_ij = NULL;
	Jacobian_s_ij = NULL;
	Jacobian_k_ij = NULL;

	MassMatrix_ij = NULL;

	mZeros_Aux = NULL;
	mId_Aux = NULL;

	Res_Stress_i = NULL;
	Res_Ext_Surf = NULL;
	Res_Time_Cont = NULL;
	Res_FSI_Cont = NULL;

	nodeReactions = NULL;

	solutionPredictor = NULL;

	normalVertex = NULL;
	stressTensor = NULL;

}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long iPoint, iElem = 0;
	unsigned short iVar, jVar, iDim, jDim, NodesElement = 0, nKindElements;

	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method
	bool fsi = config->GetFSI_Simulation();												// FSI simulation

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	double E = config->GetElasticyMod();

	nElement      = geometry->GetnElem();
	nDim          = geometry->GetnDim();
	nMarker       = geometry->GetnMarker();

	nPoint        = geometry->GetnPoint();
	nPointDomain  = geometry->GetnPointDomain();

	nKindElements = 2;

	element_container = new CElement*[nKindElements];
	node          	  = new CVariable*[nPoint];

	GradN_X = new double [nDim];
	GradN_x = new double [nDim];

	Total_CFEA			= 0.0;
	WAitken_Dyn      	= 0.0;
	WAitken_Dyn_tn1  	= 0.0;

  	SetFSI_ConvValue(0,0.0);
  	SetFSI_ConvValue(1,0.0);

	nVar = nDim;

	/*--- Define some auxiliary vectors related to the residual ---*/

	Residual = new double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
	Point_Max_Coord = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Point_Max_Coord[iVar] = new double[nDim];
		for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
	}

	/*--- Define some auxiliary vectors related to the solution ---*/

	Solution   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;

	nodeReactions = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) nodeReactions[iVar]   = 0.0;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	/*--- Check for a restart, initialize from zero otherwise ---*/

	if (!restart) {
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
			node[iPoint] = new CFEM_ElasVariable(Solution, nDim, nVar, config);
		}
	}
	else {

		/* The restart from a file needs to be implemented */

	}



	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Term ij of the Mass Matrix ---*/

	MassMatrix_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		MassMatrix_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				MassMatrix_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (constitutive contribution) ---*/

	Jacobian_c_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_c_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_c_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (stress contribution) ---*/

	Jacobian_s_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_s_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_s_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (incompressibility term) ---*/

	if (incompressible){
		Jacobian_k_ij = new double*[nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_k_ij[iVar] = new double [nVar];
				for (jVar = 0; jVar < nVar; jVar++) {
					Jacobian_k_ij[iVar][jVar] = 0.0;
				}
		}
	}
	else {
		Jacobian_k_ij = NULL;
	}

	/*--- Stress contribution to the node i ---*/
	Res_Stress_i = new double[nVar];

	/*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
	Res_Ext_Surf = new double[nVar];

	/*--- Contribution of the fluid tractions to the residual (auxiliary vector) ---*/
	if (fsi){
		Res_FSI_Cont = new double[nVar];
	}
	else {
		Res_FSI_Cont = NULL;
	}


	/*--- Time integration contribution to the residual ---*/
	if (dynamic) {
		Res_Time_Cont = new double [nVar];
	}
	else {
		Res_Time_Cont = NULL;
	}

	/*--- Matrices to impose clamped boundary conditions (TODO: Initialize them conditionally). ---*/

	mZeros_Aux = new double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mZeros_Aux[iDim] = new double[nDim];

	mId_Aux = new double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mId_Aux[iDim] = new double[nDim];

	for(iDim = 0; iDim < nDim; iDim++){
		for (jDim = 0; jDim < nDim; jDim++){
			mZeros_Aux[iDim][jDim] = 0.0;
			mId_Aux[iDim][jDim] = 0.0;
		}
		mId_Aux[iDim][iDim] = E;		// TODO: This works for clamped boundary conditions...
	}


	/*--- Initialization of matrix structures ---*/
	if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;

	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

	if (dynamic) {
		MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
		TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
		TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
	}


	/*--- Initialization of linear solver structures ---*/
	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);

	/*--- Here is where we assign the kind of each element ---*/

	if (nDim == 2){
		if (incompressible){
			element_container[EL_TRIA] = new CTRIA1(nDim, config);
			element_container[EL_QUAD] = new CQUAD4P1(nDim, config);
		}
		else{
			element_container[EL_TRIA] = new CTRIA1(nDim, config);
			element_container[EL_QUAD] = new CQUAD4(nDim, config);
		}
	}
	else if (nDim == 3){
		if (incompressible){
			element_container[EL_TETRA] = new CTETRA1(nDim, config);
			element_container[EL_HEXA] = new CHEXA8P1(nDim, config);
		}
		else{
			element_container[EL_TETRA] = new CTETRA1(nDim, config);
			element_container[EL_HEXA] = new CHEXA8(nDim, config);
		}
	}

	/*--- Initialize the auxiliary vector and matrix for the computation of the nodal Reactions ---*/

	normalVertex = new double [nDim];

	stressTensor = new double* [nDim];
	for (iVar = 0; iVar < nVar; iVar++){
		stressTensor[iVar] = new double [nDim];
	}

	/*---- Initialize the auxiliary vector for the solution predictor ---*/

	solutionPredictor = new double [nVar];

}

CFEM_ElasticitySolver::~CFEM_ElasticitySolver(void) {

	unsigned short iVar, nKindElements = 2;
	unsigned long iPoint;

	for (iPoint = 0; iPoint < nPoint; iPoint++){
		delete [] node[iPoint];
	}

	for (iVar = 0; iVar < nKindElements; iVar++){
		delete [] element_container[iVar];
	}

	for (iVar = 0; iVar < nVar; iVar++){
		delete [] Jacobian_s_ij[iVar];
		delete [] Jacobian_ij[iVar];
		delete [] Jacobian_c_ij[iVar];
		if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij[iVar];
		delete [] Point_Max_Coord[iVar];
		delete [] mZeros_Aux[iVar];
		delete [] mId_Aux[iVar];
		delete [] stressTensor[iVar];
	}

	delete [] element_container;
	delete [] node;
	delete [] Jacobian_s_ij;
	delete [] Jacobian_ij;
	delete [] Jacobian_c_ij;
	if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij;
	delete [] Res_Stress_i;
	delete [] Res_Ext_Surf;
	if (Res_Time_Cont != NULL) delete[] Res_Time_Cont;
	delete [] Solution;
	delete [] GradN_X;
	delete [] GradN_x;

	delete [] Residual;
	delete [] Residual_RMS;
	delete [] Residual_Max;
	delete [] Point_Max;
	delete [] Point_Max_Coord;

	delete [] mZeros_Aux;
	delete [] mId_Aux;

	delete [] nodeReactions;

	delete [] normalVertex;
	delete [] stressTensor;

}

void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {


	unsigned long iPoint;
	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool first_iter = (config->GetIntIter() == 0);													// Checks if it is the first iteration
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method


	/*--- Set vector entries to zero ---*/

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		LinSysAux.SetBlock_Zero(iPoint);
		LinSysRes.SetBlock_Zero(iPoint);
		LinSysSol.SetBlock_Zero(iPoint);
	}

	/*--- Set matrix entries to zero ---*/

	/*
	 * If the problem is linear, we only need one Jacobian matrix in the problem, because
	 * it is going to be constant along the calculations. Therefore, we only initialize
	 * the Jacobian matrix once, at the beginning of the simulation.
	 *
	 * We don't need first_iter, because there is only one iteration per time step in linear analysis.
	 */
	if ((initial_calc) && (linear_analysis)){
		Jacobian.SetValZero();
	}

	/*
	 * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
	 * both for linear and nonlinear analysis. Only initialized once, at the first time step.
	 *
	 * The same with the integration constants, as for now we consider the time step to be constant.
	 *
	 * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
	 */
	if ((dynamic) && (initial_calc) && (first_iter)) {
		MassMatrix.SetValZero();
		Compute_IntegrationConstants(config);
	}

	/*
	 * If the problem is nonlinear, we need to initialize the Jacobian and the stiffness matrix at least at the beginning
	 * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
	 * iteration.
	 */

	if ((nonlinear_analysis) && ((newton_raphson) || (first_iter)))	{
		Jacobian.SetValZero();
//		StiffMatrix.SetValZero();
	}

	/*
	 * Some external forces may be considered constant over the time step.
	 */
	if (first_iter)	{
		for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_SurfaceLoad_Res();
	}

}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Ks_ab;
	double *Kab = NULL;
	double *Kk_ab = NULL;
	double *Ta = NULL;
	unsigned short NelNodes, jNode;

	double checkJacobian, *checkCoord;

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {

		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
					}
				}

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);

			}

		}

	}


}

void CFEM_ElasticitySolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Ks_ab;
	double *Kab = NULL;
	double *Kk_ab = NULL;
	double *Ta = NULL;
	unsigned short NelNodes, jNode;

	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		/*--- If incompressible, we compute the Mean Dilatation term first so the volume is already computed ---*/

		if (incompressible) numerics->Compute_MeanDilatation_Term(element_container[EL_KIND]);

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);
				Ks_ab = element_container[EL_KIND]->Get_Ks_ab(iNode,jNode);
				if (incompressible) Kk_ab = element_container[EL_KIND]->Get_Kk_ab(iNode,jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					Jacobian_s_ij[iVar][iVar] = Ks_ab;
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
						if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
					}
				}

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);

			}

		}

	}

}

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Mab;
	unsigned short NelNodes, jNode;

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics->Compute_Mass_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			for (jNode = 0; jNode < NelNodes; jNode++){

				Mab = element_container[EL_KIND]->Get_Mab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					MassMatrix_ij[iVar][iVar] = Mab;
				}

				MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], MassMatrix_ij);

			}

		}

	}

}

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {


	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Ks_ab;
	double *Kab = NULL;
	double *Kk_ab = NULL;
	double *Ta = NULL;
	unsigned short NelNodes, jNode;

	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		numerics->Compute_NodalStress_Term(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

		}

	}

	for (iDim = 0; iDim < nDim; iDim++) {
		val_Coord = geometry->node[0]->GetCoord(iDim);
		val_Sol = node[0]->GetSolution(iDim) + val_Coord;
	}

}

void CFEM_ElasticitySolver::Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iDim, iStress;
	unsigned short nNodes, nStress;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	if (nDim == 2) nStress = 3;
	else if (nDim == 3) nStress = 6;

	double *Ta = NULL;

	unsigned short NelNodes;

	/*--- Restart stress to avoid adding results from previous time steps ---*/

	 for (iPoint = 0; iPoint < nPoint; iPoint++){
		for (iStress = 0; iStress < nStress; iStress++){
				node[iPoint]->SetStress_FEM(iStress, 0.0);
		}
	}

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		numerics->Compute_Averaged_NodalStress(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			/*--- This only works if the problem is nonlinear ---*/
			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysReact.AddBlock(indexNode[iNode], Res_Stress_i);

			for (iStress = 0; iStress < nStress; iStress++){
				node[indexNode[iNode]]->AddStress_FEM(iStress,
						(element_container[EL_KIND]->Get_NodalStress(iNode, iStress) /
								geometry->node[indexNode[iNode]]->GetnElem()) );
			}

		}

	}

	  double *Stress;
	  double VonMises_Stress, MaxVonMises_Stress = 0.0;
	  double Sxx,Syy,Szz,Sxy,Sxz,Syz,S1,S2;

	 /* --- For the number of nodes in the mesh ---*/
	  for (iPoint = 0; iPoint < nPoint; iPoint++) {

		  /* --- Get the stresses, added up from all the elements that connect to the node ---*/

		  Stress  = node[iPoint]->GetStress_FEM();

		  /* --- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/

		  if (geometry->GetnDim() == 2) {

			  Sxx=Stress[0];
			  Syy=Stress[1];
			  Sxy=Stress[2];

			  S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
			  S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);

			  VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);

		  }
		  else if (geometry->GetnDim() == 3) {

			  Sxx = Stress[0];
			  Syy = Stress[1];
			  Szz = Stress[3];

			  Sxy = Stress[2];
			  Sxz = Stress[4];
			  Syz = Stress[5];

			  VonMises_Stress = sqrt(0.5*(   pow(Sxx - Syy, 2.0)
											+ pow(Syy - Szz, 2.0)
											+ pow(Szz - Sxx, 2.0)
											+ 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
											));

		  }

		  node[iPoint]->SetVonMises_Stress(VonMises_Stress);

		  /*--- Compute the maximum value of the Von Mises Stress ---*/

		  MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);

	  }

  	double checkJacobian;
  	unsigned long jNode;

  	ofstream myfile;
  	myfile.open ("Reactions.txt");

  	unsigned short iMarker;
  	unsigned long iVertex;
  	double val_Reaction;

	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.

  	if (!dynamic){
  		/*--- Loop over all the markers  ---*/
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			switch (config->GetMarker_All_KindBC(iMarker)) {

				/*--- If it corresponds to a clamped boundary  ---*/

				case CLAMPED_BOUNDARY:

				myfile << "MARKER " << iMarker << ":" << endl;

					/*--- Loop over all the vertices  ---*/
					for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

					/*--- Get node index ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

					myfile << "Node " << iPoint << "." << " \t ";

					for (iDim = 0; iDim < nDim; iDim++){
						/*--- Retrieve coordinate ---*/
						val_Coord = geometry->node[iPoint]->GetCoord(iDim);
						myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
					}

					for (iVar = 0; iVar < nVar; iVar++){
						/*--- Retrieve reaction ---*/
						val_Reaction = LinSysReact.GetBlock(iPoint, iVar);
						myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
					}

					myfile << endl;
				}
			  myfile << endl;
			  break;
		}
  	}
  	else if (dynamic){

  		switch (config->GetKind_TimeIntScheme_FEA()) {
  			case (CD_EXPLICIT):
  					  cout << "NOT IMPLEMENTED YET" << endl;
  			  break;
  			case (NEWMARK_IMPLICIT):

				/*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
				if (linear_analysis){
					for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
						for (iVar = 0; iVar < nVar; iVar++){
							Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+		//a0*U(t)
										 	 a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+	//a2*U'(t)
										 	 a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
						}
						TimeRes_Aux.SetBlock(iPoint, Residual);
					}
				}
				else if (nonlinear_analysis){
					for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
						for (iVar = 0; iVar < nVar; iVar++){
							Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)  			//a0*U(t)
											 - a_dt[0]*node[iPoint]->GetSolution(iVar) 					//a0*U(t+dt)(k-1)
										 	 + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)		//a2*U'(t)
										 	 + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
						}
						TimeRes_Aux.SetBlock(iPoint, Residual);
					}
				}
				/*--- Once computed, compute M*TimeRes_Aux ---*/
				MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);

  		  		/*--- Loop over all the markers  ---*/
  				for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  					switch (config->GetMarker_All_KindBC(iMarker)) {

  						/*--- If it corresponds to a clamped boundary  ---*/

  						case CLAMPED_BOUNDARY:

  						myfile << "MARKER " << iMarker << ":" << endl;

  							/*--- Loop over all the vertices  ---*/
  							for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

  							/*--- Get node index ---*/
  							iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

  							myfile << "Node " << iPoint << "." << " \t ";

  							for (iDim = 0; iDim < nDim; iDim++){
  								/*--- Retrieve coordinate ---*/
  								val_Coord = geometry->node[iPoint]->GetCoord(iDim);
  								myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
  							}

  							/*--- Retrieve the time contribution ---*/
  							Res_Time_Cont = TimeRes.GetBlock(iPoint);

  							for (iVar = 0; iVar < nVar; iVar++){
  								/*--- Retrieve reaction ---*/
  								val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + Res_Time_Cont[iVar];
  								myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
  							}

  							myfile << endl;
  						}
  					  myfile << endl;
  					  break;
  				}


  			  break;
  			case (GA_IMPLICIT):
  		  			  cout << "NOT IMPLEMENTED YET" << endl;
  			  break;
  		  }

  	}



	myfile.close();

		#ifdef HAVE_MPI

		  /*--- Compute MaxVonMises_Stress using all the nodes ---*/

		  double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
		  MPI_Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		#endif

		  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/

	  Total_CFEA = MaxVonMises_Stress;

}

void CFEM_ElasticitySolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CFEM_ElasticitySolver::Compute_IntegrationConstants(CConfig *config) {

	double Delta_t= config->GetDelta_DynTime();
	double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();

	/*--- Integration constants for Newmark scheme ---*/

	a_dt[0]= 1 / (alpha*pow(Delta_t,2.0));
	a_dt[1]= delta / (alpha*Delta_t);
	a_dt[2]= 1 / (alpha*Delta_t);
	a_dt[3]= 1 /(2*alpha) - 1;
	a_dt[4]= delta/alpha - 1;
	a_dt[5]= (Delta_t/2) * (delta/alpha - 2);
	a_dt[6]= Delta_t * (1-delta);
	a_dt[7]= delta * Delta_t;

}


void CFEM_ElasticitySolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	unsigned short iVar, jVar;

	double tempCoord;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	unsigned short iNode, jNode;


	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		/*--- Get node index ---*/

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (nDim == 2) {
			Solution[0] = 0.0;  Solution[1] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;
		}
		else {
			Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
		}

		node[iPoint]->SetSolution(Solution);

		if (dynamic){
			node[iPoint]->SetSolution_Vel(Solution);
			node[iPoint]->SetSolution_Accel(Solution);
		}

//		for (iVar = 0; iVar < nVar; iVar++){
//			nodeReactions[iVar] = - 1.0 * LinSysRes.GetBlock(iPoint, iVar);
//		}
//
//		LinSysReact.SetBlock(iPoint,nodeReactions);

		/*--- Initialize the reaction vector ---*/
		LinSysReact.SetBlock(iPoint, Residual);


		LinSysRes.SetBlock(iPoint, Residual);

		/*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/

		/*--- Delete the columns for a particular node ---*/

		for (iVar = 0; iVar < nPoint; iVar++){
			if (iVar==iPoint) {
				Jacobian.SetBlock(iVar,iPoint,mId_Aux);
			}
			else {
				Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
			}
		}

		/*--- Delete the rows for a particular node ---*/
		for (jVar = 0; jVar < nPoint; jVar++){
			if (iPoint!=jVar) {
				Jacobian.SetBlock(iPoint,jVar,mZeros_Aux);
			}
		}

		/*--- If the problem is dynamic ---*/
		/*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/
		/*--- TODO: Do I really need to do this? ---*/

		if(dynamic){

			node[iPoint]->SetSolution_time_n(Solution);
			node[iPoint]->SetSolution_Vel_time_n(Solution);
			node[iPoint]->SetSolution_Accel_time_n(Solution);

		}

	}

}

void CFEM_ElasticitySolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		/*--- Get node index ---*/

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (nDim == 2) {
			Solution[0] = 0.0;  Solution[1] = 0.0;
		}
		else {
			Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
		}

		node[iPoint]->SetSolution(Solution);

		if (dynamic){
			node[iPoint]->SetSolution_Vel(Solution);
			node[iPoint]->SetSolution_Accel(Solution);
		}

	}

}

void CFEM_ElasticitySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
		unsigned short iMesh) {

    unsigned short iVar;
	unsigned long iPoint, total_index;

	bool first_iter = (config->GetIntIter() == 0);
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);		// Nonlinear analysis.

	double solNorm = 0.0, tempCheck[3];

	if (nonlinear_analysis){

		/*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

		/*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
		/*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
		/*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

		if (first_iter){
			Conv_Ref[0] = 1.0;											// Position for the norm of the solution
			Conv_Ref[1] = max(LinSysRes.norm(), EPS);					// Position for the norm of the residual
			Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);		// Position for the energy tolerance

			/*--- Make sure the computation runs at least 2 iterations ---*/
			Conv_Check[0] = 1.0;
			Conv_Check[1] = 1.0;
			Conv_Check[2] = 1.0;
		}
		else {
			/*--- Compute the norm of the solution vector Uk ---*/
			for (iPoint = 0; iPoint < nPoint; iPoint++){
				for (iVar = 0; iVar < nVar; iVar++){
					solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
				}
			}
			Conv_Ref[0] = max(sqrt(solNorm), EPS);							// Norm of the solution vector

			Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];					// Norm of the delta-solution vector
			Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];					// Norm of the residual
			Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];	// Position for the energy tolerance
		}

	}
	else{

			/*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

			/*---  Compute the residual Ax-f ---*/

			Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

			  /*--- Set maximum residual to zero ---*/

				for (iVar = 0; iVar < nVar; iVar++) {
					SetRes_RMS(iVar, 0.0);
					SetRes_Max(iVar, 0.0, 0);
				}

			  /*--- Compute the residual ---*/

				for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
					for (iVar = 0; iVar < nVar; iVar++) {
						total_index = iPoint*nVar+iVar;
						AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
						AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
					}
				}


			  /*--- MPI solution ---*/

			  Set_MPI_Solution(geometry, config);

			  /*--- Compute the root mean square residual ---*/

			  SetResidual_RMS(geometry, config);
	}

}

void CFEM_ElasticitySolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) {

	double a[3], b[3], AC[3], BD[3];
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	double Length_Elem = 0.0, Area_Elem = 0.0, Normal_Elem[3] = {0.0, 0.0, 0.0};
	unsigned short iDim;

	double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
	double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
	double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));

	double TotalLoad;

  bool Gradual_Load = config->GetGradual_Load();
	double CurrentTime=config->GetCurrent_DynTime();
	double ModAmpl, NonModAmpl;

  bool Ramp_Load = config->GetRamp_Load();
	double Ramp_Time = config->GetRamp_Time();

	if (Ramp_Load){
		ModAmpl=LoadDirVal*LoadDirMult*CurrentTime/Ramp_Time;
		NonModAmpl=LoadDirVal*LoadDirMult;
		TotalLoad=min(ModAmpl,NonModAmpl);
	}
	else if (Gradual_Load){
		ModAmpl=2*((1/(1+exp(-1*CurrentTime)))-0.5);
		TotalLoad=ModAmpl*LoadDirVal*LoadDirMult;
	}
	else{
		TotalLoad=LoadDirVal*LoadDirMult;
	}

	/*--- Compute the norm of the vector that was passed in the config file ---*/
	double Norm;
	if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
	if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);

	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
		if (nDim == 3) {

			Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
		    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){
		    	Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();
		    }

		}

		/*--- Compute area (3D), and length of the surfaces (2D) ---*/

		if (nDim == 2) {

			for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];

			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
			Normal_Elem[0] =   a[1];
			Normal_Elem[1] = -(a[0]);

		}

		if (nDim == 3) {

			if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){

				for (iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = Coord_1[iDim]-Coord_0[iDim];
					b[iDim] = Coord_2[iDim]-Coord_0[iDim];
				}

				double Ni=0 , Nj=0, Nk=0;

				Ni=a[1]*b[2]-a[2]*b[1];
				Nj=-a[0]*b[2]+a[2]*b[0];
				Nk=a[0]*b[1]-a[1]*b[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);

			}

			else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

				for (iDim = 0; iDim < nDim; iDim++) {
					AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
					BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
				}

				double Ni=0 , Nj=0, Nk=0;

				Ni=AC[1]*BD[2]-AC[2]*BD[1];
				Nj=-AC[0]*BD[2]+AC[2]*BD[0];
				Nk=AC[0]*BD[1]-AC[1]*BD[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);

			}
		}

      if (nDim == 2) {

        Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;

        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);

      }

      else {
    	  if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){

    		  Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

    	      node[Point_0]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_1]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_2]->Add_SurfaceLoad_Res(Residual);

    	  }
    	  else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

    		  Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

    	      node[Point_0]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_1]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_2]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_3]->Add_SurfaceLoad_Res(Residual);


    	  }

      }

	}

}

void CFEM_ElasticitySolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEM_ElasticitySolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

	unsigned long iPoint, jPoint;
	unsigned short iVar, jVar;

	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool first_iter = (config->GetIntIter() == 0);
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method
	bool fsi = config->GetFSI_Simulation();												// FSI simulation.

	double *checkCoord;

	if (!dynamic){

		for (iPoint = 0; iPoint < nPoint; iPoint++){
			/*--- Add the external contribution to the residual    ---*/
			/*--- (the terms that are constant over the time step) ---*/
			Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
		}

	}

	if (dynamic) {

		/*--- Add the mass matrix contribution to the Jacobian ---*/

		/*
		 * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
		 * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
		 * iteration, as the Jacobian is recomputed
		 *
		 * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
		 * From then on, the Jacobian is always the same matrix.
		 *
		 */

		if (((nonlinear_analysis) && ((newton_raphson) || (first_iter)))||
			((linear_analysis) && (initial_calc))) {
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
				for (jPoint = 0; jPoint < geometry->GetnPoint(); jPoint++){
					for(iVar = 0; iVar < nVar; iVar++){
						for (jVar = 0; jVar < nVar; jVar++){
							Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
						}
					}
					Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
				}
			}
		}

		/*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
		if (linear_analysis){
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				for (iVar = 0; iVar < nVar; iVar++){
					Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+		//a0*U(t)
								 	 a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+	//a2*U'(t)
								 	 a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
				}
				TimeRes_Aux.SetBlock(iPoint, Residual);
			}
		}
		else if (nonlinear_analysis){
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				for (iVar = 0; iVar < nVar; iVar++){
					Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)  			//a0*U(t)
									 - a_dt[0]*node[iPoint]->GetSolution(iVar) 					//a0*U(t+dt)(k-1)
								 	 + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)		//a2*U'(t)
								 	 + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
				}
				TimeRes_Aux.SetBlock(iPoint, Residual);
			}
		}
		/*--- Once computed, compute M*TimeRes_Aux ---*/
		MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
		/*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			/*--- Dynamic contribution ---*/
			Res_Time_Cont = TimeRes.GetBlock(iPoint);
			LinSysRes.AddBlock(iPoint, Res_Time_Cont);
			/*--- External surface load contribution ---*/
			Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
			checkCoord = geometry->node[iPoint]->GetCoord();
			/*--- Add FSI contribution ---*/
			if (fsi) {
				/*--- TODO: It may be worthy restricting the flow traction to the boundary elements... ---*/
				Res_FSI_Cont = node[iPoint]->Get_FlowTraction();
				LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
			}
		}
	}



}

void CFEM_ElasticitySolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

    unsigned short iVar;
	unsigned long iPoint, total_index;

	bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);		// Geometrically linear problems
	bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.

	unsigned short iNode, jNode, jVar;

	/*--- Update solution ---*/

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

		for (iVar = 0; iVar < nVar; iVar++) {

			/*--- Displacements component of the solution ---*/

			/*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/

			if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

			if (nonlinear)	node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

		}

	}

	if (dynamic){

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

			for (iVar = 0; iVar < nVar; iVar++) {

				/*--- Acceleration component of the solution ---*/
				/*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

				Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
										node[iPoint]->GetSolution_time_n(iVar)) -
							   a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
							   a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);

			}

			/*--- Set the acceleration in the node structure ---*/

			node[iPoint]->SetSolution_Accel(Solution);

			for (iVar = 0; iVar < nVar; iVar++) {

				/*--- Velocity component of the solution ---*/
				/*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

				Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
							   a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
							   a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);

			}

			/*--- Set the velocity in the node structure ---*/

			node[iPoint]->SetSolution_Vel(Solution);

		}

	}


}

void CFEM_ElasticitySolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config){

	unsigned long IterLinSol;

	CSysSolve femSystem;
	IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

}



void CFEM_ElasticitySolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry,
											CGeometry **flow_geometry, CConfig *fea_config,
											CConfig *flow_config, CNumerics *fea_numerics) {



	unsigned short nMarkerFSIint, nMarkerFEA, nMarkerFlow;		// Number of markers on FSI problem, FEA and Flow side
	unsigned short iMarkerFSIint, iMarkerFEA, iMarkerFlow;		// Variables for iteration over markers
	unsigned short markFEA, markFlow;

	unsigned long nVertexFEA, nVertexFlow;						// Number of vertices on FEA and Flow side
	unsigned long iVertex, iPoint;								// Variables for iteration over vertices and nodes


	/*--- TODO: We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/
	/*--- This may be improved ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++){
		node[iPoint]->Clear_FlowTraction();
	}

	unsigned short iDim, jDim;

	// Check the kind of fluid problem
	bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
	bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
								(flow_config->GetKind_Solver() == RANS) );

	unsigned long nodeVertex, donorVertex;
	double *normalsVertex;

	double *tn_f;
	tn_f 				= new double [nVar];			// Fluid traction

  	/*--- Redimensionalize the pressure ---*/

	double *Velocity_ND, *Velocity_Real;
	double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
	double factorForces;

    Velocity_Real = flow_config->GetVelocity_FreeStream();
    Density_Real = flow_config->GetDensity_FreeStream();

    Velocity_ND = flow_config->GetVelocity_FreeStreamND();
    Density_ND = flow_config->GetDensity_FreeStreamND();

	Velocity2_Real = 0.0;
	Velocity2_ND = 0.0;
    for (iDim = 0; iDim < nDim; iDim++){
    	Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    	Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }

    factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

	/*--- Apply a ramp to the transfer of the fluid loads ---*/

	double ModAmpl;
	double CurrentTime=fea_config->GetCurrent_DynTime();
	double Static_Time=fea_config->GetStatic_Time();

	bool Ramp_Load = fea_config->GetRamp_Load();
	double Ramp_Time = fea_config->GetRamp_Time();

	if (CurrentTime <= Static_Time){ ModAmpl=0.0; }
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Ramp_Time)) &&
			(Ramp_Load)){
		ModAmpl=(CurrentTime-Static_Time)/Ramp_Time;
		ModAmpl=max(ModAmpl,0.0);
		ModAmpl=min(ModAmpl,1.0);
	}
	else{ ModAmpl=1.0; }

	// Parameters for the calculations
	// Pn: Pressure
	// Pinf: Pressure_infinite
	// div_vel: Velocity divergence
	// Dij: Dirac delta
	double Pn = 0.0, Pinf = 0.0, div_vel = 0.0, Dij = 0.0;
	double Viscosity = 0.0, Density = 0.0;
	double **Grad_PrimVar;
	double Tau[3][3];

	/*--- Number of markers in the FSI interface ---*/
	nMarkerFSIint = (fea_config->GetMarker_n_FSIinterface())/2;

	nMarkerFEA  = fea_geometry[MESH_0]->GetnMarker();		// Retrieve total number of markers on FEA side
	nMarkerFlow = flow_geometry[MESH_0]->GetnMarker();		// Retrieve total number of markers on Fluid side

	/*--- Loop over all the markers on the interface ---*/

	for (iMarkerFSIint=0; iMarkerFSIint < nMarkerFSIint; iMarkerFSIint++){

		/*--- Identification of the markers ---*/

		/*--- Current structural marker ---*/
		for (iMarkerFEA=0; iMarkerFEA < nMarkerFEA; iMarkerFEA++){
			if ( fea_config->GetMarker_All_FSIinterface(iMarkerFEA) == (iMarkerFSIint+1)){
				markFEA=iMarkerFEA;
			}
		}

		/*--- Current fluid marker ---*/
		for (iMarkerFlow=0; iMarkerFlow < nMarkerFlow; iMarkerFlow++){
			if (flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == (iMarkerFSIint+1)){
				markFlow=iMarkerFlow;
			}
		}

		nVertexFEA = fea_geometry[MESH_0]->GetnVertex(markFEA);		// Retrieve total number of vertices on FEA marker
		nVertexFlow = flow_geometry[MESH_0]->GetnVertex(markFlow);  // Retrieve total number of vertices on Fluid marker

		/*--- Loop over the nodes in the structural mesh, calculate the tf vector (unitary) ---*/
		/*--- Here, we are looping over the fluid, and we find the pointer to the structure (donorVertex) ---*/
		for (iVertex=0; iVertex < nVertexFlow; iVertex++){

			// Node from the flow mesh
			nodeVertex=flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetNode();

			// Normals at the vertex: these normals go inside the fluid domain.
			normalsVertex = flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetNormal();

			// Corresponding node on the structural mesh
			donorVertex = flow_geometry[MESH_0]->vertex[markFlow][iVertex]->GetDonorPoint();

			// Retrieve the values of pressure, viscosity and density
			if (incompressible){

				Pn=flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetPressureInc();
				Pinf=flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetLaminarViscosityInc();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetDensityInc();

				}
			}
			else if (compressible){

				Pn=flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetPressure();
				Pinf=flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetLaminarViscosity();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[nodeVertex]->GetDensity();

				}
			}

			// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
			for (iDim = 0; iDim < nDim; iDim++) {
				tn_f[iDim] = -(Pn-Pinf)*normalsVertex[iDim];
			}

			// Calculate tn in the fluid nodes for the viscous term

			if (viscous_flow){

				// Divergence of the velocity
				div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
				if (incompressible) div_vel = 0.0;

				for (iDim = 0; iDim < nDim; iDim++) {

					for (jDim = 0 ; jDim < nDim; jDim++) {
						// Dirac delta
						Dij = 0.0; if (iDim == jDim) Dij = 1.0;

						// Viscous stress
						Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
								TWO3*Viscosity*div_vel*Dij;

						// Viscous component in the tn vector --> Units of force (non-dimensional).
						tn_f[iDim] += Tau[iDim][jDim]*normalsVertex[jDim];
					}
				}
			}

			// Rescale tn to SI units and apply time-dependent coefficient (static structure, ramp load, full load)

			for (iDim = 0; iDim < nDim; iDim++) {
				Residual[iDim] = tn_f[iDim]*factorForces*ModAmpl;
			}

			/*--- Set the Flow traction ---*/
			//node[donorVertex]->Set_FlowTraction(Residual);
			/*--- Add to the Flow traction (to add values to corners...) ---*/
			node[donorVertex]->Add_FlowTraction(Residual);
		}

	}


}

void CFEM_ElasticitySolver::SetFEA_Load_Int(CSolver ***flow_solution, CGeometry **fea_geometry,
											CGeometry **flow_geometry, CConfig *fea_config,
											CConfig *flow_config, CNumerics *fea_numerics){ }

void CFEM_ElasticitySolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                            				CConfig *fea_config, CSolver ***fea_solution){

    unsigned short predOrder = fea_config->GetPredictorOrder();
	double Delta_t = fea_config->GetDelta_DynTime();
    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    double *solDisp, *solVel, *solVel_tn, *valPred;
    double *DisplacementDonor, *SolutionDonor;

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();


    for (iPoint=0; iPoint < nPoint; iPoint++){
    	if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    	else if (predOrder==1) {

    		solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
    		solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
    		valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

    		for (iDim=0; iDim<nDim; iDim++){
    			valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
    		}

    	}
    	else if (predOrder==2) {

    		solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
    		solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
    		solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
    		valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

    		for (iDim=0; iDim<nDim; iDim++){
    			valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
    		}

    	}
    	else {
    		cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
    		fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    	}
    }


}

void CFEM_ElasticitySolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
        				  CSolver ***fea_solution, unsigned long iFSIIter){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
    double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
    double delta_deltaU[3] = {0.0, 0.0, 0.0};
    double numAitk, denAitk, WAitken;
	double CurrentTime=fea_config->GetCurrent_DynTime();
	double Static_Time=fea_config->GetStatic_Time();
	double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn;

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    WAitken=fea_config->GetAitkenStatRelax();

	numAitk = 0.0;
	denAitk = 0.0;

	ofstream historyFile_FSI;
	bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
	if (writeHistFSI){
		char cstrFSI[200];
		string filenameHistFSI = fea_config->GetConv_FileName_FSI();
		strcpy (cstrFSI, filenameHistFSI.data());
		historyFile_FSI.open (cstrFSI, std::ios_base::app);
	}


	/*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/

	if (CurrentTime > Static_Time) {

		if (iFSIIter == 0){

			WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
			WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();

			WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);

			/*--- Temporal fix, only for now ---*/
			WAitkDyn = max(WAitkDyn, 0.1);

			SetWAitken_Dyn(WAitkDyn);
			if (writeHistFSI){
				historyFile_FSI << " " << endl ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
				historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
			}

		}
		else{

			for (iPoint=0; iPoint<nPoint; iPoint++){

				dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
				dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
				dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
				dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();

				for (iDim=0; iDim < nDim; iDim++){

					/*--- Compute the deltaU and deltaU_n+1 ---*/
					deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
					deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];

					/*--- Compute the difference ---*/
					delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];

					/*--- Add numerator and denominator ---*/
					numAitk += deltaU[iDim] * delta_deltaU[iDim];
					denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];

				}

			}

				WAitkDyn = GetWAitken_Dyn();

			if (denAitk > 1E-8){
				WAitkDyn = - 1.0 * WAitkDyn * numAitk / denAitk ;
			}

				WAitkDyn = max(WAitkDyn, 0.1);
				WAitkDyn = min(WAitkDyn, 1.0);

				SetWAitken_Dyn(WAitkDyn);

				if (writeHistFSI){
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
					historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
				}

		}

	}

	if (writeHistFSI){historyFile_FSI.close();}

}

void CFEM_ElasticitySolver::SetAitken_Relaxation(CGeometry **fea_geometry,
        				  CConfig *fea_config, CSolver ***fea_solution){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    unsigned short RelaxMethod_FSI;
    double *dispPred, *dispCalc;
    double WAitken;
	double CurrentTime=fea_config->GetCurrent_DynTime();
	double Static_Time=fea_config->GetStatic_Time();

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();

	/*--- Only when there is movement it makes sense to update the solutions... ---*/

	if (CurrentTime > Static_Time) {

		if (RelaxMethod_FSI == NO_RELAXATION){
			WAitken = 1.0;
		}
		else if (RelaxMethod_FSI == FIXED_PARAMETER){
			WAitken = fea_config->GetAitkenStatRelax();
		}
		else if (RelaxMethod_FSI == AITKEN_DYNAMIC){
			WAitken = GetWAitken_Dyn();
		}
		else {
			WAitken = 1.0;
			cout << "No relaxation parameter used. " << endl;
		}


		for (iPoint=0; iPoint<nPoint; iPoint++){

			/*--- Retrieve pointers to the predicted and calculated solutions ---*/
			dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
			dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();

			/*--- Set predicted solution as the old predicted solution ---*/
			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();

			/*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);

			/*--- Apply the Aitken relaxation ---*/
			for (iDim=0; iDim < nDim; iDim++){
				dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
			}

		}

	}

}

void CFEM_ElasticitySolver::Update_StructSolution(CGeometry **fea_geometry,
        				  CConfig *fea_config, CSolver ***fea_solution){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    double *valSolutionPred, *valSolution;

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    for (iPoint=0; iPoint<nPoint; iPoint++){

    	valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

		fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);

    }

}
