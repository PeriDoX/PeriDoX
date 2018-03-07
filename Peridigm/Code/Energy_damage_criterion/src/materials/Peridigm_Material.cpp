/*! \file Peridigm_Material.cpp */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "Peridigm_Material.hpp"
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
//#include "material_utilities.h"
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::SerialMatrix& jacobian,
                                           PeridigmNS::Material::JacobianType jacobianType) const
{
  // Compute a finite-difference Jacobian using either FORWARD_DIFFERENCE or CENTRAL_DIFFERENCE
  computeFiniteDifferenceJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, CENTRAL_DIFFERENCE, jacobianType);
}

void PeridigmNS::Material::computeFiniteDifferenceJacobian(const double dt,
                                                           const int numOwnedPoints,
                                                           const int* ownedIDs,
                                                           const int* neighborhoodList,
                                                           PeridigmNS::DataManager& dataManager,
                                                           PeridigmNS::SerialMatrix& jacobian,
                                                           FiniteDifferenceScheme finiteDifferenceScheme,
                                                           PeridigmNS::Material::JacobianType jacobianType) const
{
  // The Jacobian is of the form:
  //
  // dF_0x/dx_0  dF_0x/dy_0  dF_0x/dz_0  dF_0x/dx_1  dF_0x/dy_1  dF_0x/dz_1  ...  dF_0x/dx_n  dF_0x/dy_n  dF_0x/dz_n  
  // dF_0y/dx_0  dF_0y/dy_0  dF_0y/dz_0  dF_0y/dx_1  dF_0y/dy_1  dF_0y/dz_1  ...  dF_0y/dx_n  dF_0y/dy_n  dF_0y/dz_n  
  // dF_0z/dx_0  dF_0z/dy_0  dF_0z/dz_0  dF_0z/dx_1  dF_0z/dy_1  dF_0z/dz_1  ...  dF_0z/dx_n  dF_0z/dy_n  dF_0z/dz_n  
  // dF_1x/dx_0  dF_1x/dy_0  dF_1x/dz_0  dF_1x/dx_1  dF_1x/dy_1  dF_1x/dz_1  ...  dF_1x/dx_n  dF_1x/dy_n  dF_1x/dz_n  
  // dF_1y/dx_0  dF_1y/dy_0  dF_1y/dz_0  dF_1y/dx_1  dF_1y/dy_1  dF_1y/dz_1  ...  dF_1y/dx_n  dF_1y/dy_n  dF_1y/dz_n  
  // dF_1z/dx_0  dF_1z/dy_0  dF_1z/dz_0  dF_1z/dx_1  dF_1z/dy_1  dF_1z/dz_1  ...  dF_1z/dx_n  dF_1z/dy_n  dF_1z/dz_n  
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  // dF_nx/dx_0  dF_nx/dy_0  dF_nx/dz_0  dF_nx/dx_1  dF_nx/dy_1  dF_nx/dz_1  ...  dF_nx/dx_n  dF_nx/dy_n  dF_nx/dz_n  
  // dF_ny/dx_0  dF_ny/dy_0  dF_ny/dz_0  dF_ny/dx_1  dF_ny/dy_1  dF_ny/dz_1  ...  dF_ny/dx_n  dF_ny/dy_n  dF_ny/dz_n  
  // dF_nz/dx_0  dF_nz/dy_0  dF_nz/dz_0  dF_nz/dx_1  dF_nz/dy_1  dF_nz/dz_1  ...  dF_nz/dx_n  dF_nz/dy_n  dF_nz/dz_n  

  // Each entry is computed by finite difference:
  //
  // Forward difference:
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon
  //
  // Central difference:
  // dF_0x/dx_0 = ( F_0x(positive perturbed x_0) - F_0x(negative perturbed x_0) ) / ( 2.0*epsilon )

  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_finiteDifferenceProbeLength == DBL_MAX, "**** Finite-difference Jacobian requires that the \"Finite Difference Probe Length\" parameter be set.\n");
  double epsilon = m_finiteDifferenceProbeLength;

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId = fieldManager.getFieldId("Volume");
  int coordinatesFId = fieldManager.getFieldId("Coordinates");
  int velocityFId = fieldManager.getFieldId("Velocity");
  int forceDensityFId = fieldManager.getFieldId("Force_Density");

  // Loop over all points.
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors+1);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numNeighbors+1); 
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same fields and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data in the constitutiveData array.
    double *volume, *y, *v, *force;
    tempDataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    tempDataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(velocityFId, PeridigmField::STEP_NP1)->ExtractView(&v);
    tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1)->ExtractView(&force);

    // Create a temporary vector for storing force
    Teuchos::RCP<Epetra_Vector> forceVector = tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1);
    Teuchos::RCP<Epetra_Vector> tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
    double* tempForce;
    tempForceVector->ExtractView(&tempForce);

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < 3*(numNeighbors+1))
      scratchMatrix.Resize(3*(numNeighbors+1));

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(3*(numNeighbors+1));
    for(int i=0 ; i<numNeighbors+1 ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<3 ; ++j)
        globalIndices[3*i+j] = 3*globalID+j;
    }

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE){
      // Compute and store the unperturbed force.
      computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      for(int i=0 ; i<forceVector->MyLength() ; ++i)
        tempForce[i] = force[i];
    }

    // Perturb one dof in the neighborhood at a time and compute the force.
    // The point itself plus each of its neighbors must be perturbed.
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;

      for(int dof=0 ; dof<3 ; ++dof){

        // Perturb a dof and compute the forces.
        double oldY = y[3*perturbID+dof];
        double oldV = v[3*perturbID+dof];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // Compute and store the negatively perturbed force.
          y[3*perturbID+dof] -= epsilon;
          v[3*perturbID+dof] -= epsilon/dt;
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[3*perturbID+dof] = oldY;
          v[3*perturbID+dof] = oldV;
          for(int i=0 ; i<forceVector->MyLength() ; ++i)
            tempForce[i] = force[i];
        }


        // Compute the purturbed force
        y[3*perturbID+dof] += epsilon;
        v[3*perturbID+dof] += epsilon/dt;
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[3*perturbID+dof] = oldY;
        v[3*perturbID+dof] = oldV;

        for(int i=0 ; i<numNeighbors+1 ; ++i){
          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;

          for(int d=0 ; d<3 ; ++d){
            double value = ( force[3*forceID+d] - tempForce[3*forceID+d] ) / epsilon;
            if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(3*forceID+d, 3*perturbID+dof) = value;
          }
        }
      }
    }

    // Convert force density to force
    // \todo Create utility function for this in ScratchMatrix
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        scratchMatrix(row, col) *= volume[row/3];
      }
    }

    // Check for NaNs
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(scratchMatrix(row, col)), "**** NaN detected in finite-difference Jacobian.\n");
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
      jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL) {
      jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    }
    else // unknown jacobian type
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Unknown Jacobian Type\n");
  }
}

double PeridigmNS::Material::calculateBulkModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) + 
    static_cast<int>(shearModulusDefined) + 
    static_cast<int>(youngsModulusDefined) + 
    static_cast<int>(poissonsRatioDefined);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(bulkModulusDefined)
    computedValue = bulkModulus;
  else if(youngsModulusDefined && shearModulusDefined)
    computedValue = (youngsModulus * shearModulus) / (3.0*(3.0*shearModulus - youngsModulus));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (3.0*(1.0 - 2.0*poissonsRatio));
  else if(shearModulusDefined && poissonsRatioDefined)
    computedValue = (2.0*shearModulus*(1.0 + poissonsRatio)) / (3.0*(1.0 - 2.0*poissonsRatio));
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}


double PeridigmNS::Material::calculateShearModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) + 
    static_cast<int>(shearModulusDefined) + 
    static_cast<int>(youngsModulusDefined) + 
    static_cast<int>(poissonsRatioDefined);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(shearModulusDefined)
    computedValue = shearModulus;
  else if(bulkModulusDefined && youngsModulusDefined)
    computedValue = (3.0*bulkModulus*youngsModulus) / (9.0*bulkModulus - youngsModulus);
  else if(bulkModulusDefined & poissonsRatioDefined)
    computedValue = (3.0*bulkModulus*(1.0 - 2.0*poissonsRatio)) / (2.0*(1.0 + poissonsRatio));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (2.0*(1.0 + poissonsRatio));
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}


void
PeridigmNS::Material::computeDilatation(
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        const double BM,
        const double SM,
        PeridigmNS::DataManager& dataManager) const {
            
          
    //////////////////////////////////////////////////////////////////////////////////////    
    // Routine developed for Peridigm by DLR.
    // Dilatation has to be updated before the damage routine. Therefore, this routine has been added
    // to the model evaluator. The dilatation is needed to calculate the energy at both bonds (12 and 21) 
    // for the correct time step, prior to the calculation of the internal forces. The bonds with critical
    // energy are deleted and the internal loads calculation will be done.
    //////////////////////////////////////////////////////////////////////////////////////    
    // Questions?
    // Christian Willberg  christian.willberg@dlr.de
    // Martin Raedel       martin.raedel@dlr.de
    //////////////////////////////////////////////////////////////////////////////////////   
     
    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    //double m_horizon = params.get<double>("Horizon");
    //////////////////////////////////////////////////////////////////////////////////////
    // ALL IDS will be initialized in the material routine "Peridigm_ElasticMaterial.cpp"
    // carefull if changed from Peridigm_ElasticMaterial.cpp to something else
    //////////////////////////////////////////////////////////////////////////////////////    
    
    int m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
    int m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    int m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    int m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
    int m_dilatationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
	int m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
    int m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
    int m_damageModelFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Damage_Model_Data");
   
    ///////////////////////////////////////////////////////////////////////////////
    PeridigmNS::InfluenceFunction::functionPointer OMEGA;
    OMEGA = PeridigmNS::InfluenceFunction::self().getInfluenceFunction();
    double *x, *y, *bondDamage, *horizon;
    double *Volume, *weightedVolume, *damageModel, *theta;
    ///////////////////////////////////////////////////////////////////////////////
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&Volume);
	dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);   
    dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&theta);
	dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    
    dataManager.getData(m_damageModelFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageModel);

    
	*(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
    double m_alpha = 0.0; // to be adapted for temperature
    double *deltaTemperature = NULL; 
	const double *xOwned = x;
	const double *yOwned = y;
	const double *deltaT = deltaTemperature;
	const double *m = weightedVolume;
	const double *v = Volume;
	double cellVolume;
	
	const int *neighPtr = neighborhoodList;
	
	
	for(int p=0; p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, horizon++, deltaT++, m++, theta++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;
		*theta = 0.0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &x[3*localId];
			const double *YP = &y[3*localId];
			double X_dx = XP[0]-X[0];
			double X_dy = XP[1]-X[1];
			double X_dz = XP[2]-X[2];
			double zetaSquared = X_dx*X_dx+X_dy*X_dy+X_dz*X_dz;
			double Y_dx = YP[0]-Y[0];
			double Y_dy = YP[1]-Y[1];
			double Y_dz = YP[2]-Y[2];
			double dY = Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz;
			double d = sqrt(zetaSquared);
			double e = sqrt(dY);
			e -= d;
			//if(deltaTemperature)
			//  e -= thermalExpansionCoefficient*(*deltaT)*d;
			double omega = OMEGA(d,*horizon);
			*theta += 3.0*omega*(1.0-*bondDamage)*d*e*cellVolume/(*m);
		}
        // variable is needed to determine the accurate energy in the energy based damage criterion
		damageModel[3*p] = *theta;
		damageModel[3*p+1] = BM / (*m);
		damageModel[3*p+2] = SM / (*m) ;
		*theta = 0.0;

	}
}		

