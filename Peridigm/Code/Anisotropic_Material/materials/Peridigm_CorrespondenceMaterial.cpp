/*! \file Peridigm_CorrespondenceMaterial.cpp */

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

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include "elastic_correspondence.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
using namespace std;

PeridigmNS::CorrespondenceMaterial::CorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), m_hourglassCoefficient(0.0),
        m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_shapeTensorInverseFieldId(-1),
    m_leftStretchTensorFieldId(-1),
    m_rotationTensorFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1),
    m_cauchyStressFieldId(-1), 
    m_unrotatedRateOfDeformationFieldId(-1),
    m_partialStressFieldId(-1),
    m_hourglassStiffId(-1)
    
{
     
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  
  m_stabilizationType = 3;

  m_plane = false;

  //************************************
  // wie komme ich an den Namen??
  //************************************
  //if (params.isParameter("Linear Elastic Correspondence")){
    nonLin = false;
  
    bool m_planeStrain = false, m_planeStress = false;
    if (params.isParameter("Plane Strain"))
        m_planeStrain = params.get<bool>("Plane Strain");
        
    if (params.isParameter("Plane Stress"))
        m_planeStress = params.get<bool>("Plane Stress");

    if (m_planeStrain==true){
        m_plane=true;
        
    }
    if (m_planeStress==true){
        m_plane=true;
       
    }
  m_tension = true;
  if (params.isParameter("Tension pressure separation for damage model")){
      m_tension = params.get<bool>("Tension pressure separation for damage model");
      
  }
  m_applyAutomaticDifferentiationJacobian = true;
  if(params.isParameter("Apply Automatic Differentiation Jacobian"))
    m_applyAutomaticDifferentiationJacobian = params.get<bool>("Apply Automatic Differentiation Jacobian");
  if (params.isParameter("Stabilizaton Type")){
      
    if (params.get<string>("Stabilizaton Type")=="Bond Based"){
        m_stabilizationType = 1;
        m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    }
    if (params.get<string>("Stabilizaton Type")=="State Based"){
        m_stabilizationType = 2;
        m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    }
    //if (params.get<string>("Stabilizaton Type")=="Sub Horizon"){
    //    m_stabilizationType = 4;
    //    // works only for linear elastic correspondence first
    //    // based on: Shubhankar Roy Chowdhury, Pranesh Roy, Debasish Roy and J N Reddy, 
    //    // "A simple alteration of the peridynamics correspondence principle to eliminate zero-energy deformation"
    //}
    if (params.get<string>("Stabilizaton Type")=="Global Stiffness"){
        
        m_stabilizationType = 3;
        m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");

        double C11 = 0.0, C12 = 0.0, C13 = 0.0, C14 = 0.0, C15 = 0.0, C16 = 0.0;
        double            C22 = 0.0, C23 = 0.0, C24 = 0.0, C25 = 0.0, C26 = 0.0;
        double                       C33 = 0.0, C34 = 0.0, C35 = 0.0, C36 = 0.0;
        double                                  C44 = 0.0, C45 = 0.0, C46 = 0.0;
        double                                             C55 = 0.0, C56 = 0.0;
        double                                                        C66 = 0.0;
       
        bool iso = false;
        if (params.isParameter("Material Symmetry")){
            if (params.get<string>("Material Symmetry")=="Isotropic"){
                C11 = params.get<double>("C11");
                C44 = params.get<double>("C44");
                C55 = params.get<double>("C44");
                C66 = params.get<double>("C44");
                C12 = C11 - 2*C55;
                C13 = C12;
                C14 = 0.0;
                C15 = 0.0;
                C16 = 0.0;
                C22 = params.get<double>("C11");
                C33 = params.get<double>("C11");
                C23 = C12;
                C24 = 0.0;
                C25 = 0.0;
                C26 = 0.0;
                C34 = 0.0;
                C35 = 0.0;
                C36 = 0.0;
                C45 = 0.0;
                C46 = 0.0;
                C56 = 0.0; 
                iso  = true;
            }
            if (params.get<string>("Material Symmetry")=="Anisotropic"){
                C11 = params.get<double>("C11");
                C12 = params.get<double>("C12");
                C13 = params.get<double>("C13");
                C14 = params.get<double>("C14");
                C15 = params.get<double>("C15");
                C16 = params.get<double>("C16");
                C22 = params.get<double>("C22");
                C23 = params.get<double>("C23");
                C24 = params.get<double>("C24");
                C25 = params.get<double>("C25");
                C26 = params.get<double>("C26");
                C33 = params.get<double>("C33");
                C34 = params.get<double>("C34");
                C35 = params.get<double>("C35");
                C36 = params.get<double>("C36");
                C44 = params.get<double>("C44");
                C45 = params.get<double>("C45");
                C46 = params.get<double>("C46");
                C55 = params.get<double>("C55");
                C56 = params.get<double>("C56");
                C66 = params.get<double>("C66");
                
            }
        
        }
        else{
            m_bulkModulus = calculateBulkModulus(params);
            m_shearModulus = calculateShearModulus(params);
        
            C11 = (4*m_shearModulus*(3*m_bulkModulus + m_shearModulus))/(3*m_bulkModulus + 4*m_shearModulus);
            C44 = m_shearModulus;
            C55 = m_shearModulus;
            C66 = m_shearModulus;
            C12 = (2*(3*m_bulkModulus - 2*m_shearModulus)*m_shearModulus)/(3*m_bulkModulus + 4*m_shearModulus);
            C13 = C12;
            C14 = 0.0;
            C15 = 0.0;
            C16 = 0.0;
            C22 = C11;
            C33 = C11;
            C23 = C12;
            C24 = 0.0;
            C25 = 0.0;
            C26 = 0.0;
            C34 = 0.0;
            C35 = 0.0;
            C36 = 0.0;
            C45 = 0.0;
            C46 = 0.0;
            C56 = 0.0;
        }
        // Equation (8) Dipasquale, D., Sarego, G., Zaccariotto, M., Galvanetto, U., A discussion on failure criteria
          // for ordinary state-based Peridynamics, Engineering Fracture Mechanics (2017), doi: https://doi.org/10.1016/
          // j.engfracmech.2017.10.011
        if (m_planeStrain==true)m_plane=true;
        if (m_planeStress==true)m_plane=true;
        // have to be done after rotation if angles exist
        if (m_plane==false){
         C[0][0] = C11;C[0][1] = C12;C[0][2]= C13; C[0][3] = C14; C[0][4] = C15; C[0][5]= C16;
         C[1][0] = C12;C[1][1] = C22;C[1][2]= C23; C[1][3] = C24; C[1][4] = C25; C[1][5]= C26;
         C[2][0] = C13;C[2][1] = C23;C[2][2]= C33; C[2][3] = C34; C[2][4] = C35; C[2][5]= C36;
         C[3][0] = C14;C[3][1] = C24;C[3][2]= C34; C[3][3] = C44; C[3][4] = C45; C[3][5]= C46;
         C[4][0] = C15;C[4][1] = C25;C[4][2]= C35; C[4][3] = C45; C[4][4] = C55; C[4][5]= C56;
         C[5][0] = C16;C[5][1] = C26;C[5][2]= C36; C[5][3] = C46; C[5][4] = C56; C[5][5]= C66;
        }
        // tbd in future
        if (m_planeStress==true && iso == true){
            //only transversal isotropic in the moment --> definition of iso missing
         C[0][0] = C11-C13*C13/C22;C[0][1] = C12-C13*C23/C22;C[0][2] = 0.0; C[0][3] = 0.0; C[0][4] = 0.0; C[0][5] = 0.0;
         C[1][0] = C12-C13*C23/C22;C[1][1] = C22-C13*C23/C22;C[1][2] = 0.0; C[1][3] = 0.0; C[1][4] = 0.0; C[1][5] = 0.0;
         C[2][0] = 0.0;            C[2][1] = 0.0;            C[2][2] = 0.0; C[2][3] = 0.0; C[2][4] = 0.0; C[2][5] = 0.0;
         C[3][0] = 0.0;            C[3][1] = 0.0;            C[3][2] = 0.0; C[3][3] = 0.0; C[3][4] = 0.0; C[3][5] = 0.0;
         C[4][0] = 0.0;            C[4][1] = 0.0;            C[4][2] = 0.0; C[4][3] = 0.0; C[4][4] = 0.0; C[4][5] = 0.0;
         C[5][0] = 0.0;            C[5][1] = 0.0;            C[5][2] = 0.0; C[5][3] = 0.0; C[5][4] = 0.0; C[5][5] = C66;
        }
        // not correct for plane stress!!!
        if (m_plane==true && iso == false){
         C[0][0] = C11;C[0][1] = C12;C[0][2] = 0.0;C[0][3] = 0.0;C[0][4] = 0.0;C[0][5] = C16;
         C[1][0] = C12;C[1][1] = C22;C[1][2] = 0.0;C[1][3] = 0.0;C[1][4] = 0.0;C[1][5] = C26;
         C[2][0] = 0.0;C[2][1] = 0.0;C[2][2] = 0.0;C[2][3] = 0.0;C[2][4] = 0.0;C[2][5] = 0.0;
         C[3][0] = 0.0;C[3][1] = 0.0;C[3][2] = 0.0;C[3][3] = 0.0;C[3][4] = 0.0;C[3][5] = 0.0;
         C[4][0] = 0.0;C[4][1] = 0.0;C[4][2] = 0.0;C[4][3] = 0.0;C[4][4] = 0.0;C[4][5] = 0.0;
         C[5][0] = C16;C[5][1] = C26;C[5][2] = 0.0;C[5][3] = 0.0;C[5][4] = 0.0;C[5][5] = C66;
        
        }
        
    }
  }
  nonLin = false;
  lin = true;

  if (params.isParameter("Elastic Correspondence")){
     
      nonLin = true;
      lin = false;

    }




  //TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_modelAnglesId                     = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_hourglassForceDensityFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Hourglass_Force_Density");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
  m_leftStretchTensorFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_shapeTensorInverseFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_piolaStressTimesInvShapeTensorXId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorX");
  m_piolaStressTimesInvShapeTensorYId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorY");
  m_piolaStressTimesInvShapeTensorZId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorZ");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  m_detachedNodesFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
  m_hourglassStiffId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Hourglass_Stiffness");
  
  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_hourglassForceDensityFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorXId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorYId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorZId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_detachedNodesFieldId);
  m_fieldIds.push_back(m_hourglassStiffId);
  m_fieldIds.push_back(m_modelAnglesId);

  }

PeridigmNS::CorrespondenceMaterial::~CorrespondenceMaterial()
{
}

void
PeridigmNS::CorrespondenceMaterial::initialize(const double dt,
                                               const int numOwnedPoints,
                                               const int* ownedIDs,
                                               const int* neighborhoodList,
                                               PeridigmNS::DataManager& dataManager)
{
  
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  
  
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  
  dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);


  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  double *detachedNodes;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  //Initialize the left stretch and rotation tenor to the identity matrix
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorNP1, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorNP1, numOwnedPoints);

  //Initialize the inverse of the shape tensor and the deformation gradient
  double *volume;
  double *horizon;
  double *modelCoordinates;
  double *coordinates;
  double *coordinatesNP1;
  double *shapeTensorInverse;
  double *deformationGradient;
  double *bondDamage;


  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&coordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

 int shapeTensorReturnCode = 0;
  shapeTensorReturnCode = 
            CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                   horizon,
                                                                                   modelCoordinates,
                                                                                   coordinates,
                                                                                   coordinatesNP1,
                                                                                   shapeTensorInverse,
                                                                                   deformationGradient,
                                                                                   bondDamage,
                                                                                   neighborhoodList,
                                                                                   numOwnedPoints,
                                                                                   m_plane,
                                                                                   detachedNodes);

  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::initialize() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

}

void
PeridigmNS::CorrespondenceMaterial::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces and partial stress

  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  double *horizon, *volume, *modelCoordinates, *coordinates, *coordinatesNP1, *shapeTensorInverse, *deformationGradient, *bondDamage, *pointAngles, *detachedNodes;
  //double *deformationGradientNonInc;

  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_modelAnglesId,           PeridigmField::STEP_NONE)->ExtractView(&pointAngles);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&coordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);

          
          
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation

   // vier szenarien
   int shapeTensorReturnCode = 0;
  // if (lin == true){
    
    
    shapeTensorReturnCode = 
        CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                horizon,
                                                                                modelCoordinates,
                                                                                coordinates,
                                                                                coordinatesNP1,
                                                                                shapeTensorInverse,
                                                                                deformationGradient,
                                                                                bondDamage,
                                                                                neighborhoodList,
                                                                                numOwnedPoints,
                                                                                m_plane,
                                                                                detachedNodes);

    
   //}

  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial::computeForce() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

  
  double *velocities, *leftStretchTensorN, *leftStretchTensorNP1, *rotationTensorN, *rotationTensorNP1, *unrotatedRateOfDeformation;
  if (nonLin==true){  
      dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
      dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
      dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
      dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
      dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);
      dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
      // Compute left stretch tensor, rotation tensor, and unrotated rate-of-deformation.
      // Performs a polar decomposition via Flanagan & Taylor (1987) algorithm.
      //"A non-ordinary state-based peridynamic method to model solid material deformation and fracture"
      // Thomas L. Warren, Stewart A. Silling, Abe Askari, Olaf Weckner, Michael A. Epton, Jifeng Xu
      
      int rotationTensorReturnCode = 0;

      rotationTensorReturnCode = CORRESPONDENCE::computeUnrotatedRateOfDeformationAndRotationTensor(volume,
                                                                                                   horizon,
                                                                                                   modelCoordinates, 
                                                                                                   velocities, 
                                                                                                   deformationGradient,
                                                                                                   shapeTensorInverse,
                                                                                                   leftStretchTensorN,
                                                                                                   rotationTensorN,
                                                                                                   leftStretchTensorNP1,
                                                                                                   rotationTensorNP1,
                                                                                                   unrotatedRateOfDeformation,
                                                                                                   neighborhoodList, 
                                                                                                   numOwnedPoints, 
                                                                                                   dt,
                                                                                                   bondDamage,
                                                                                                   m_plane,
                                                                                                   detachedNodes);

      string rotationTensorErrorMessage =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
      rotationTensorErrorMessage +=
        "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
      string rotationTensorErrorMessage2 =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to invert deformation gradient tensor.\n";
        
      string rotationTensorErrorMessage3 =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to invert temp.\n";

      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 1, rotationTensorErrorMessage);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 2, rotationTensorErrorMessage2);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode == 3, rotationTensorErrorMessage3);
                                    
     }  
    
  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  // computeCauchyStress() typically uses the following fields which are accessed via the DataManager:
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  
  // multiple Cauchy stresses will be provided over the datamanager --> Peridigm_ElasticLinearCorrespondence

                
  computeCauchyStress(dt, numOwnedPoints, dataManager);

  // rotate back to the Eulerian frame
  double *unrotatedCauchyStressNP1, *cauchyStressNP1;


  if (nonLin == true) {
      dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
      dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);
      
      CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                         unrotatedCauchyStressNP1,
                                         cauchyStressNP1,
                                         numOwnedPoints);
                                         // Cauchy stress is now updated and in the rotated state. 
  }
  else
  {
    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);
    
  }

   
  // Proceed with conversion to Piola-Kirchoff and force-vector states.
//--------------------------------------------------------------------------
 double *forceDensity;
 dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);
 dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
 double *partialStress;
 dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);
 dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
 
 double *tempStressX, *tempStressY, *tempStressZ;
 dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->ExtractView(&tempStressX);
 dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->ExtractView(&tempStressY);
 dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->ExtractView(&tempStressZ);
 dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
 dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
 dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);

 double *hourglassStiff;

 dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->ExtractView(&hourglassStiff);



 CORRESPONDENCE::computeForcesAndStresses(
                                       numOwnedPoints,
                                       neighborhoodList,
                                       volume,
                                       horizon,
                                       modelCoordinates,
                                       coordinatesNP1,
                                       deformationGradient,
                                       cauchyStressNP1,
                                       shapeTensorInverse,
                                       bondDamage,
                                       C,
                                       pointAngles,
                                       forceDensity,
                                       partialStress,
                                       tempStressX,
                                       tempStressY,
                                       tempStressZ,
                                       hourglassStiff,
                                       m_hourglassCoefficient,
                                       m_stabilizationType,
                                       m_plane,
                                       m_tension,
                                       detachedNodes);
                                          
 //     std::cout<<numOwnedPoints<< " "<< *(deformationGradient)<<" "<<*(partialStress)<<std::endl;
    

  if (m_incremental == false){
      if (m_stabilizationType == 1){
      CORRESPONDENCE::computeHourglassForce(volume,
                                            horizon,
                                            modelCoordinates,
                                            coordinates,
                                            deformationGradient,
                                            forceDensity,
                                            neighborhoodList,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_hourglassCoefficient,
                                            bondDamage
                                            );
      }
      
      if (m_stabilizationType == 2){
      
      CORRESPONDENCE::computeCorrespondenceStabilityForce(volume,
                                            horizon,
                                            modelCoordinates,
                                            coordinates,
                                            deformationGradient,
                                            forceDensity,
                                            neighborhoodList,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_hourglassCoefficient,
                                            bondDamage
                                            );
      }
  }
}

void
PeridigmNS::CorrespondenceMaterial::computeJacobian(const double dt,
                                             const int numOwnedPoints,
                                             const int* ownedIDs,
                                             const int* neighborhoodList,
                                             PeridigmNS::DataManager& dataManager,
                                             PeridigmNS::SerialMatrix& jacobian,
                                             PeridigmNS::Material::JacobianType jacobianType) const
{

  if(m_applyAutomaticDifferentiationJacobian){
    // Compute the Jacobian via automatic differentiation
    computeAutomaticDifferentiationJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, jacobianType);  
  }
  else{
  //  // Call the base class function, which computes the Jacobian by finite difference
    computeJacobianFiniteDifference(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, CENTRAL_DIFFERENCE, jacobianType);
  }
}


void
PeridigmNS::CorrespondenceMaterial::computeAutomaticDifferentiationJacobian(const double dt,
                                                                     const int numOwnedPoints,
                                                                     const int* ownedIDs,
                                                                     const int* neighborhoodList,
                                                                     PeridigmNS::DataManager& dataManager,
                                                                     PeridigmNS::SerialMatrix& jacobian,
                                                                     PeridigmNS::Material::JacobianType jacobianType) const
{
  // Compute contributions to the tangent matrix on an element-by-element basis

  // To reduce memory re-allocation, use static variable to store Fad types for
  // current coordinates (independent variables).

  static vector<Sacado::Fad::DFad<double> > coordinatesNP1_AD;
  static vector<Sacado::Fad::DFad<double> > coordinates_AD;

  static Sacado::Fad::DFad<double> C_AD[6][6];

  // Loop over all points.
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    int numEntries = numNeighbors+1;
    int numDof = 3*numEntries;
    vector<int> tempMyGlobalIDs(numEntries);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numEntries); 
    tempNeighborhoodList[0] = numNeighbors;

    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same field specs and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(tempNumOwnedPoints);
    tempOwnedIDs[0] = 0;

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < numDof)
      scratchMatrix.Resize(numDof);

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof);
    for(int i=0 ; i<numEntries ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<3 ; ++j)
        globalIndices[3*i+j] = 3*globalID+j;
    }

    // Extract pointers to the underlying data in the constitutive data array.

    double *modelCoordinates, *coordinates, *volume, *coordinatesNP1, *angles, *detachedNodes;
    tempDataManager.getData(m_modelAnglesId,           PeridigmField::STEP_NONE)->ExtractView(&angles);
    tempDataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    tempDataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
    tempDataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&coordinates);
    tempDataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
    double *horizon, *bondDamage;
    tempDataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    tempDataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
 
    // Create arrays of Fad objects for the current coordinates, dilatation, and force density
    // Modify the existing vector of Fad objects for the current coordinates
    if((int)coordinatesNP1_AD.size() < numDof) {
      coordinates_AD.resize(numDof);
      coordinatesNP1_AD.resize(numDof);
    }
    for(int i=0 ; i<numDof ; ++i){
      coordinates_AD[i].diff(i, numDof);
      coordinates_AD[i].val() = coordinates[i];
      coordinatesNP1_AD[i].diff(i, numDof);
      coordinatesNP1_AD[i].val() = coordinatesNP1[i];
    }
    
    vector<Sacado::Fad::DFad<double> > shapeTensorInverse_AD;
    vector<Sacado::Fad::DFad<double> > deformationGradient_AD;
    vector<Sacado::Fad::DFad<double> > cauchyStress_AD;
    vector<Sacado::Fad::DFad<double> > cauchyStressNP1_AD;
    vector<Sacado::Fad::DFad<double> > partialStress_AD;
    vector<Sacado::Fad::DFad<double> > tempStressX_AD;
    vector<Sacado::Fad::DFad<double> > tempStressY_AD;
    vector<Sacado::Fad::DFad<double> > tempStressZ_AD;
    vector<Sacado::Fad::DFad<double> > hourglassStiff_AD;
    
    partialStress_AD.resize(numDof*numDof);
    shapeTensorInverse_AD.resize(numDof*numDof);
    deformationGradient_AD.resize(numDof*numDof);
    cauchyStress_AD.resize(numDof*numDof);
    cauchyStressNP1_AD.resize(numDof*numDof);
    tempStressX_AD.resize(numDof*numDof);
    tempStressY_AD.resize(numDof*numDof);
    tempStressZ_AD.resize(numDof*numDof);
    hourglassStiff_AD.resize(numDof*numDof);

    for(int i=0 ; i<6 ; ++i){
        for(int j=0 ; j<6 ; ++j){
            C_AD[i][j].val() = C[i][j];
    }}

    int shapeTensorReturnCode = 0;
      shapeTensorReturnCode = 
                CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                       horizon,
                                                                                       modelCoordinates,
                                                                                       &coordinates_AD[0],
                                                                                       &coordinatesNP1_AD[0],
                                                                                       &shapeTensorInverse_AD[0],
                                                                                       &deformationGradient_AD[0],
                                                                                       bondDamage,
                                                                                       &tempNeighborhoodList[0],
                                                                                       tempNumOwnedPoints,
                                                                                       m_plane,
                                                                                       detachedNodes);
    double *cauchyStress;
    double *cauchyStressNP1;

    tempDataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStress);
    tempDataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNP1);
    CORRESPONDENCE::updateElasticCauchyStressSmallDef(&deformationGradient_AD[0], 
                                          &cauchyStress_AD[0],
                                          &cauchyStressNP1_AD[0],
                                          tempNumOwnedPoints,
                                          &C_AD[0],
                                          angles,
                                          m_type,
                                          dt,
                                          m_incremental);
    
    // define the sacadoFAD force vector
    vector<Sacado::Fad::DFad<double> > force_AD(numDof);
    //std::cout<<m_stabilizationType<<std::endl;
    CORRESPONDENCE::computeForcesAndStresses(
                                        tempNumOwnedPoints,
                                        &tempNeighborhoodList[0],
                                        volume,
                                        horizon,
                                        modelCoordinates,
                                        &coordinatesNP1_AD[0],
                                        &deformationGradient_AD[0],
                                        &cauchyStressNP1_AD[0],
                                        &shapeTensorInverse_AD[0],
                                        bondDamage,
                                        &C_AD[0],
                                        angles,
                                        &force_AD[0],
                                        &partialStress_AD[0],
                                        &tempStressX_AD[0],
                                        &tempStressY_AD[0],
                                        &tempStressZ_AD[0],
                                        &hourglassStiff_AD[0],
                                        m_hourglassCoefficient,
                                        m_stabilizationType,
                                        m_plane,
                                        m_tension,
                                        detachedNodes);

    // Load derivative values into scratch matrix
    // Multiply by volume along the way to convert force density to force
    double value;
    for(int row=0 ; row<numDof ; ++row){
      for(int col=0 ; col<numDof ; ++col){
	value = force_AD[row].dx(col) ; //--> I think this must be it, because forces are already provided
    //value = force_AD[row].dx(col) * volume[row/3]; // given by peridigm org
	TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(value), "**** NaN detected in correspondence::computeAutomaticDifferentiationJacobian().\n");
        scratchMatrix(row, col) = value;
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

void
PeridigmNS::CorrespondenceMaterial::computeJacobianFiniteDifference(const double dt,
                                                           const int numOwnedPoints,
                                                           const int* ownedIDs,
                                                           const int* neighborhoodList,
                                                           PeridigmNS::DataManager& dataManager,
                                                           PeridigmNS::SerialMatrix& jacobian,
                                                           FiniteDifferenceScheme finiteDifferenceScheme,
                                                           PeridigmNS::Material::JacobianType jacobianType) const
{

  
  // The mechanics Jacobian is of the form:
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

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  bool solveForDisplacement = dofManager.displacementTreatedAsUnknown();
  //bool solveForTemperature = dofManager.temperatureTreatedAsUnknown();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();
  int numDisplacementDof = dofManager.numberOfDisplacementDegreesOfFreedom();
  //int numTemperatureDof = dofManager.numberOfTemperatureDegreesOfFreedom();
  int displacementDofOffset = dofManager.displacementDofOffset();
  //int temperatureDofOffset = dofManager.temperatureDofOffset();

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId(-1), coordinatesFId(-1), velocityFId(-1), forceDensityFId(-1); //, temperatureFId(-1), fluxDivergenceFId(-1);
  volumeFId = fieldManager.getFieldId("Volume");
  if (solveForDisplacement) {
    coordinatesFId = fieldManager.getFieldId("Coordinates");
    velocityFId = fieldManager.getFieldId("Velocity");
    forceDensityFId = fieldManager.getFieldId("Force_Density");
  }
  //if (solveForTemperature) {
  //  temperatureFId = fieldManager.getFieldId("Temperature");
  //  fluxDivergenceFId = fieldManager.getFieldId("Flux_Divergence");
  //}

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

    // Extract pointers to the underlying data.
    double *volume, *y, *v, *force;//, *temperature, *fluxDivergence;
    tempDataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    if (solveForDisplacement) {
      tempDataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
      tempDataManager.getData(velocityFId, PeridigmField::STEP_NP1)->ExtractView(&v);
      tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1)->ExtractView(&force);
    }
    //if (solveForTemperature) {
    //  tempDataManager.getData(temperatureFId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
    //  tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1)->ExtractView(&fluxDivergence);
    //}

    // Create a temporary vector for storing force and/or flux divergence.
    Teuchos::RCP<Epetra_Vector> forceVector, tempForceVector, fluxDivergenceVector, tempFluxDivergenceVector;
    double *tempForce, *tempFluxDivergence;
    if (solveForDisplacement) {
      forceVector = tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1);
      tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
      tempForceVector->ExtractView(&tempForce);
    }
    //if (solveForTemperature) {
    //  fluxDivergenceVector = tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1);
    //  tempFluxDivergenceVector = Teuchos::rcp(new Epetra_Vector(*fluxDivergenceVector));
    //  tempFluxDivergenceVector->ExtractView(&tempFluxDivergence);
    //}

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < numDof*(numNeighbors+1))
      scratchMatrix.Resize(numDof*(numNeighbors+1));

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof*(numNeighbors+1));
    for(int i=0 ; i<numNeighbors+1 ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<numDof ; ++j){
        globalIndices[numDof*i+j] = numDof*globalID+j;
      }
    }

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE){
      if (solveForDisplacement) {
        // Compute and store the unperturbed force.
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        
        
        
        for(int i=0 ; i<forceVector->MyLength() ; ++i)
          tempForce[i] = force[i];
      }
      //if (solveForTemperature) {
      //  // Compute and store the unperturbed flux divergence.
      //  computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      //  for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
      //    tempFluxDivergence[i] = fluxDivergence[i];
      //}
    }

    // Perturb one dof in the neighborhood at a time and compute the force and/or flux divergence.
    // The point itself plus each of its neighbors must be perturbed.
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;

      // Displacement degrees of freedom
      for(int dof=0 ; dof<numDisplacementDof ; ++dof){

        // Perturb a dof and compute the forces.
        double oldY = y[numDof*perturbID+dof];
        double oldV = v[numDof*perturbID+dof];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // Compute and store the negatively perturbed force.
          y[numDof*perturbID+dof] -= epsilon;
          v[numDof*perturbID+dof] -= epsilon/dt;
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[numDof*perturbID+dof] = oldY;
          v[numDof*perturbID+dof] = oldV;
          for(int i=0 ; i<forceVector->MyLength() ; ++i)
            tempForce[i] = force[i];
        }

        // Compute the purturbed force.
        y[numDof*perturbID+dof] += epsilon;
        v[numDof*perturbID+dof] += epsilon/dt;
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[numDof*perturbID+dof] = oldY;
        v[numDof*perturbID+dof] = oldV;

        for(int i=0 ; i<numNeighbors+1 ; ++i){
          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;

          for(int d=0 ; d<numDof ; ++d){
            double value = ( force[numDof*forceID+d] - tempForce[numDof*forceID+d] ) / epsilon;
            if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(numDof*forceID + displacementDofOffset + d, numDof*perturbID + displacementDofOffset + dof) = value;
          }
        }
      }

      // Temperature degrees of freedom
    //if(solveForTemperature){
    //
    //  // Perturb a temperature value and compute the flux divergence.
    //  double oldTemperature = temperature[perturbID];
    //
    //  if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
    //    // Compute and store the negatively perturbed flux divergence.
    //    temperature[perturbID] -= epsilon;
    //    computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
    //    temperature[perturbID] = oldTemperature;
    //    for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
    //      tempFluxDivergence[i] = fluxDivergence[i];
    //  }
    //
    //// Compute the purturbed flux divergence.
    //  temperature[perturbID] += epsilon;
    //  computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
    //  temperature[perturbID] = oldTemperature;
    //
    //  for(int i=0 ; i<numNeighbors+1 ; ++i){
    //    int fluxDivergenceID;
    //    if(i < numNeighbors)
    //      fluxDivergenceID = tempNeighborhoodList[i+1];
    //    else
    //      fluxDivergenceID = 0;
    //
    //    double value = ( fluxDivergence[fluxDivergenceID] - tempFluxDivergence[fluxDivergenceID] ) / epsilon;
    //    if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
    //      value *= 0.5;
    //    scratchMatrix(numDof*fluxDivergenceID + temperatureDofOffset, numDof*perturbID + temperatureDofOffset) = value;
    //}
    //}
    }

    // Multiply by nodal volume
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        scratchMatrix(row, col) *= volume[row/numDof];
      }
    }

    // Check for NaNs
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(scratchMatrix(row, col)), "**** NaN detected in finite-difference Jacobian.\n");
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
