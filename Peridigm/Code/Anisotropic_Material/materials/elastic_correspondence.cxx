//! \file elastic_correspondence.cxx

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

#include "elastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h> // for sqrt

//#include <Teuchos_Assert.hpp>
//#include <Epetra_SerialComm.h>


using namespace std;
namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticCauchyStress
(
const ScalarT* deltaTemperatureN,
const ScalarT* deltaTemperatureNP1,
const ScalarT* unrotatedRateOfDeformation,
const ScalarT* unrotatedCauchyStressN,
ScalarT* unrotatedCauchyStressNP1,
const int numPoints,
const double bulkMod,
const double shearMod,
const double alpha,
const double dt
)
{
  const ScalarT* deltaTempN = deltaTemperatureN;
  const ScalarT* deltaTempNP1 = deltaTemperatureNP1;

  // Hooke's law
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  for(int iID=0 ; iID<numPoints ; ++iID,
      rateOfDef+=9, sigmaN+=9, sigmaNP1+=9){

      //strainInc = dt * rateOfDef
      for(int i = 0; i < 9; i++){
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //thermal strains
      if(deltaTemperatureN && deltaTemperatureNP1){
        double thermalStrainN = alpha*deltaTemperatureN[iID];
        double thermalStrainNP1 = alpha*deltaTemperatureNP1[iID];
        dilatationInc -= 3.0*(thermalStrainNP1 - thermalStrainN);
      }

      //update stress
      for(int i = 0; i < 9; i++){
          *(sigmaNP1+i) = *(sigmaN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(sigmaNP1) += bulkMod*dilatationInc;
      *(sigmaNP1+4) += bulkMod*dilatationInc;
      *(sigmaNP1+8) += bulkMod*dilatationInc;

  }
}

template<typename ScalarT>
void updateElasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* unrotatedCauchyStressN, 
    ScalarT* unrotatedCauchyStressNP1, 
    ScalarT* vonMisesStress,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;

  ScalarT* vmStress = vonMisesStress;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  ScalarT deviatoricStressNP1[9];
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;

  for(int iID=0 ; iID<numPoints ; ++iID, 
      rateOfDef+=9, sigmaN+=9, sigmaNP1+=9, ++vmStress){

    //strainInc = dt * rateOfDef
    for(int i = 0; i < 9; i++){
        strainInc[i] = *(rateOfDef+i)*dt;
        deviatoricStrainInc[i] = strainInc[i];
    }

    //dilatation
    dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

    //deviatoric strain
    deviatoricStrainInc[0] -= dilatationInc/3.0;
    deviatoricStrainInc[4] -= dilatationInc/3.0;
    deviatoricStrainInc[8] -= dilatationInc/3.0;

    //update stress
    for(int i = 0; i < 9; i++){
      *(sigmaNP1+i) = *(sigmaN+i) + deviatoricStrainInc[i]*2.0*shearMod;
    }
    *(sigmaNP1) += bulkMod*dilatationInc;
    *(sigmaNP1+4) += bulkMod*dilatationInc;
    *(sigmaNP1+8) += bulkMod*dilatationInc;

    sphericalStressNP1 = (*(sigmaNP1) + *(sigmaNP1+4) + *(sigmaNP1+8))/3.0;

    // Compute the ``trial'' von Mises stress
    for(int i = 0; i < 9; i++){
      deviatoricStressNP1[i] = *(sigmaNP1+i);
    }
    deviatoricStressNP1[0] -= sphericalStressNP1;
    deviatoricStressNP1[4] -= sphericalStressNP1;
    deviatoricStressNP1[8] -= sphericalStressNP1;

    // Compute \sigma_ij * \sigma_ij
    tempScalar = 0.0;
    for(int j = 0; j < 3; j++){
      for(int i = 0; i < 3; i++){
        tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
      }
    }

    *vmStress = sqrt(3.0/2.0*tempScalar);
  }
}

template<typename ScalarT>
void updateElasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* unrotatedCauchyStressN, 
    ScalarT* unrotatedCauchyStressNP1, 
    ScalarT* vonMisesStress,
    const double* flyingPointFlag,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;

  ScalarT* vmStress = vonMisesStress;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  ScalarT deviatoricStressNP1[9];
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;

  for(int iID=0 ; iID<numPoints ; ++iID, 
      rateOfDef+=9, sigmaN+=9, sigmaNP1+=9, ++vmStress, ++flyingPointFlg){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      //strainInc = dt * rateOfDef
      for(int i = 0; i < 9; i++){
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //update stress
      for(int i = 0; i < 9; i++){
          *(sigmaNP1+i) = *(sigmaN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(sigmaNP1) += bulkMod*dilatationInc;
      *(sigmaNP1+4) += bulkMod*dilatationInc;
      *(sigmaNP1+8) += bulkMod*dilatationInc;

      sphericalStressNP1 = (*(sigmaNP1) + *(sigmaNP1+4) + *(sigmaNP1+8))/3.0;

      // Compute the ``trial'' von Mises stress
      for(int i = 0; i < 9; i++){
          deviatoricStressNP1[i] = *(sigmaNP1+i);
      }
      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
              tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
      }

      *vmStress = sqrt(3.0/2.0*tempScalar);
    }
  }

}

template<typename ScalarT>
void updateBondLevelElasticCauchyStress
(
    const ScalarT* bondLevelUnrotatedRateOfDeformationXX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZZ, 
    const ScalarT* bondLevelUnrotatedCauchyStressXXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZZN, 
    ScalarT* bondLevelUnrotatedCauchyStressXXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZZNP1, 
    ScalarT* bondLevelVonMisesStress,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  const ScalarT* rateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  const ScalarT* rateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  const ScalarT* rateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  const ScalarT* rateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  const ScalarT* rateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  const ScalarT* rateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  const ScalarT* rateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  const ScalarT* rateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;
  const ScalarT* sigmaXXN = bondLevelUnrotatedCauchyStressXXN;
  const ScalarT* sigmaXYN = bondLevelUnrotatedCauchyStressXYN;
  const ScalarT* sigmaXZN = bondLevelUnrotatedCauchyStressXZN;
  const ScalarT* sigmaYXN = bondLevelUnrotatedCauchyStressYXN;
  const ScalarT* sigmaYYN = bondLevelUnrotatedCauchyStressYYN;
  const ScalarT* sigmaYZN = bondLevelUnrotatedCauchyStressYZN;
  const ScalarT* sigmaZXN = bondLevelUnrotatedCauchyStressZXN;
  const ScalarT* sigmaZYN = bondLevelUnrotatedCauchyStressZYN;
  const ScalarT* sigmaZZN = bondLevelUnrotatedCauchyStressZZN;
  ScalarT* sigmaXXNP1 = bondLevelUnrotatedCauchyStressXXNP1;
  ScalarT* sigmaXYNP1 = bondLevelUnrotatedCauchyStressXYNP1;
  ScalarT* sigmaXZNP1 = bondLevelUnrotatedCauchyStressXZNP1;
  ScalarT* sigmaYXNP1 = bondLevelUnrotatedCauchyStressYXNP1;
  ScalarT* sigmaYYNP1 = bondLevelUnrotatedCauchyStressYYNP1;
  ScalarT* sigmaYZNP1 = bondLevelUnrotatedCauchyStressYZNP1;
  ScalarT* sigmaZXNP1 = bondLevelUnrotatedCauchyStressZXNP1;
  ScalarT* sigmaZYNP1 = bondLevelUnrotatedCauchyStressZYNP1;
  ScalarT* sigmaZZNP1 = bondLevelUnrotatedCauchyStressZZNP1;
  
  ScalarT* vmStress = bondLevelVonMisesStress;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  ScalarT deviatoricStressNP1[9];
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

    // All is bond level.
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
        rateOfDefXX++, rateOfDefXY++, rateOfDefXZ++, 
        rateOfDefYX++, rateOfDefYY++, rateOfDefYZ++, 
        rateOfDefZX++, rateOfDefZY++, rateOfDefZZ++,
        sigmaXXN++, sigmaXYN++, sigmaXZN++, 
        sigmaYXN++, sigmaYYN++, sigmaYZN++, 
        sigmaZXN++, sigmaZYN++, sigmaZZN++,
        sigmaXXNP1++, sigmaXYNP1++, sigmaXZNP1++, 
        sigmaYXNP1++, sigmaYYNP1++, sigmaYZNP1++, 
        sigmaZXNP1++, sigmaZYNP1++, sigmaZZNP1++,
        vmStress++){

      //strainInc = dt * rateOfDef
      strainInc[0] = *rateOfDefXX*dt; strainInc[1] = *rateOfDefXY*dt; strainInc[2] = *rateOfDefXZ*dt;
      strainInc[3] = *rateOfDefYX*dt; strainInc[4] = *rateOfDefYY*dt; strainInc[5] = *rateOfDefYZ*dt;
      strainInc[6] = *rateOfDefZX*dt; strainInc[7] = *rateOfDefZY*dt; strainInc[8] = *rateOfDefZZ*dt;

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      for (int i = 0; i < 9; i++) {
        deviatoricStrainInc[i] = strainInc[i];
      }
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //update stress
      *sigmaXXNP1 = *sigmaXXN + deviatoricStrainInc[0]*2.0*shearMod;
      *sigmaXYNP1 = *sigmaXYN + deviatoricStrainInc[1]*2.0*shearMod;
      *sigmaXZNP1 = *sigmaXZN + deviatoricStrainInc[2]*2.0*shearMod;
      *sigmaYXNP1 = *sigmaYXN + deviatoricStrainInc[3]*2.0*shearMod;
      *sigmaYYNP1 = *sigmaYYN + deviatoricStrainInc[4]*2.0*shearMod;
      *sigmaYZNP1 = *sigmaYZN + deviatoricStrainInc[5]*2.0*shearMod;
      *sigmaZXNP1 = *sigmaZXN + deviatoricStrainInc[6]*2.0*shearMod;
      *sigmaZYNP1 = *sigmaZYN + deviatoricStrainInc[7]*2.0*shearMod;
      *sigmaZZNP1 = *sigmaZZN + deviatoricStrainInc[8]*2.0*shearMod;

      *sigmaXXNP1 += bulkMod*dilatationInc;
      *sigmaYYNP1 += bulkMod*dilatationInc;
      *sigmaZZNP1 += bulkMod*dilatationInc;

      // compute von mises stress for failure analysis
      sphericalStressNP1 = (*sigmaXXNP1 + *sigmaYYNP1 + *sigmaZZNP1)/3.0;

      // Compute the ``trial'' von Mises stress
      deviatoricStressNP1[0] = *sigmaXXNP1; deviatoricStressNP1[1] = *sigmaXYNP1; deviatoricStressNP1[2] = *sigmaXZNP1;
      deviatoricStressNP1[3] = *sigmaYXNP1; deviatoricStressNP1[4] = *sigmaYYNP1; deviatoricStressNP1[5] = *sigmaYZNP1;
      deviatoricStressNP1[6] = *sigmaZXNP1; deviatoricStressNP1[7] = *sigmaZYNP1; deviatoricStressNP1[8] = *sigmaZZNP1;

      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
        for(int i = 0; i < 3; i++){
          tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
        }
      }

      *vmStress = sqrt(3.0/2.0*tempScalar);
    }
  }
}

template<typename ScalarT>
void updateBondLevelElasticCauchyStress
(
    const ScalarT* bondLevelUnrotatedRateOfDeformationXX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZZ, 
    const ScalarT* bondLevelUnrotatedCauchyStressXXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZZN, 
    ScalarT* bondLevelUnrotatedCauchyStressXXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZZNP1, 
    ScalarT* bondLevelVonMisesStress,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  const ScalarT* rateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  const ScalarT* rateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  const ScalarT* rateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  const ScalarT* rateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  const ScalarT* rateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  const ScalarT* rateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  const ScalarT* rateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  const ScalarT* rateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;
  const ScalarT* sigmaXXN = bondLevelUnrotatedCauchyStressXXN;
  const ScalarT* sigmaXYN = bondLevelUnrotatedCauchyStressXYN;
  const ScalarT* sigmaXZN = bondLevelUnrotatedCauchyStressXZN;
  const ScalarT* sigmaYXN = bondLevelUnrotatedCauchyStressYXN;
  const ScalarT* sigmaYYN = bondLevelUnrotatedCauchyStressYYN;
  const ScalarT* sigmaYZN = bondLevelUnrotatedCauchyStressYZN;
  const ScalarT* sigmaZXN = bondLevelUnrotatedCauchyStressZXN;
  const ScalarT* sigmaZYN = bondLevelUnrotatedCauchyStressZYN;
  const ScalarT* sigmaZZN = bondLevelUnrotatedCauchyStressZZN;
  ScalarT* sigmaXXNP1 = bondLevelUnrotatedCauchyStressXXNP1;
  ScalarT* sigmaXYNP1 = bondLevelUnrotatedCauchyStressXYNP1;
  ScalarT* sigmaXZNP1 = bondLevelUnrotatedCauchyStressXZNP1;
  ScalarT* sigmaYXNP1 = bondLevelUnrotatedCauchyStressYXNP1;
  ScalarT* sigmaYYNP1 = bondLevelUnrotatedCauchyStressYYNP1;
  ScalarT* sigmaYZNP1 = bondLevelUnrotatedCauchyStressYZNP1;
  ScalarT* sigmaZXNP1 = bondLevelUnrotatedCauchyStressZXNP1;
  ScalarT* sigmaZYNP1 = bondLevelUnrotatedCauchyStressZYNP1;
  ScalarT* sigmaZZNP1 = bondLevelUnrotatedCauchyStressZZNP1;
  
  ScalarT* vmStress = bondLevelVonMisesStress;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  ScalarT deviatoricStressNP1[9];
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          rateOfDefXX++, rateOfDefXY++, rateOfDefXZ++, 
          rateOfDefYX++, rateOfDefYY++, rateOfDefYZ++, 
          rateOfDefZX++, rateOfDefZY++, rateOfDefZZ++,
          sigmaXXN++, sigmaXYN++, sigmaXZN++, 
          sigmaYXN++, sigmaYYN++, sigmaYZN++, 
          sigmaZXN++, sigmaZYN++, sigmaZZN++,
          sigmaXXNP1++, sigmaXYNP1++, sigmaXZNP1++, 
          sigmaYXNP1++, sigmaYYNP1++, sigmaYZNP1++, 
          sigmaZXNP1++, sigmaZYNP1++, sigmaZZNP1++,
          vmStress++){

        //strainInc = dt * rateOfDef
        strainInc[0] = *rateOfDefXX*dt; strainInc[1] = *rateOfDefXY*dt; strainInc[2] = *rateOfDefXZ*dt;
        strainInc[3] = *rateOfDefYX*dt; strainInc[4] = *rateOfDefYY*dt; strainInc[5] = *rateOfDefYZ*dt;
        strainInc[6] = *rateOfDefZX*dt; strainInc[7] = *rateOfDefZY*dt; strainInc[8] = *rateOfDefZZ*dt;

        //dilatation
        dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

        //deviatoric strain
        for(int i = 0; i < 9; i++){
          deviatoricStrainInc[i] = strainInc[i];
        }
        deviatoricStrainInc[0] -= dilatationInc/3.0;
        deviatoricStrainInc[4] -= dilatationInc/3.0;
        deviatoricStrainInc[8] -= dilatationInc/3.0;

        //update stress
        *sigmaXXNP1 = *sigmaXXN + deviatoricStrainInc[0]*2.0*shearMod;
        *sigmaXYNP1 = *sigmaXYN + deviatoricStrainInc[1]*2.0*shearMod;
        *sigmaXZNP1 = *sigmaXZN + deviatoricStrainInc[2]*2.0*shearMod;
        *sigmaYXNP1 = *sigmaYXN + deviatoricStrainInc[3]*2.0*shearMod;
        *sigmaYYNP1 = *sigmaYYN + deviatoricStrainInc[4]*2.0*shearMod;
        *sigmaYZNP1 = *sigmaYZN + deviatoricStrainInc[5]*2.0*shearMod;
        *sigmaZXNP1 = *sigmaZXN + deviatoricStrainInc[6]*2.0*shearMod;
        *sigmaZYNP1 = *sigmaZYN + deviatoricStrainInc[7]*2.0*shearMod;
        *sigmaZZNP1 = *sigmaZZN + deviatoricStrainInc[8]*2.0*shearMod;

        *sigmaXXNP1 += bulkMod*dilatationInc;
        *sigmaYYNP1 += bulkMod*dilatationInc;
        *sigmaZZNP1 += bulkMod*dilatationInc;

        // compute von mises stress for failure analysis
        sphericalStressNP1 = (*sigmaXXNP1 + *sigmaYYNP1 + *sigmaZZNP1)/3.0;

        // Compute the ``trial'' von Mises stress
        deviatoricStressNP1[0] = *sigmaXXNP1; deviatoricStressNP1[1] = *sigmaXYNP1; deviatoricStressNP1[2] = *sigmaXZNP1;
        deviatoricStressNP1[3] = *sigmaYXNP1; deviatoricStressNP1[4] = *sigmaYYNP1; deviatoricStressNP1[5] = *sigmaYZNP1;
        deviatoricStressNP1[6] = *sigmaZXNP1; deviatoricStressNP1[7] = *sigmaZYNP1; deviatoricStressNP1[8] = *sigmaZZNP1;

        deviatoricStressNP1[0] -= sphericalStressNP1;
        deviatoricStressNP1[4] -= sphericalStressNP1;
        deviatoricStressNP1[8] -= sphericalStressNP1;

        // Compute \sigma_ij * \sigma_ij
        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
        }

        *vmStress = sqrt(3.0/2.0*tempScalar);
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      rateOfDefXX += numNeighbors; rateOfDefXY += numNeighbors; rateOfDefXZ += numNeighbors; 
      rateOfDefYX += numNeighbors; rateOfDefYY += numNeighbors; rateOfDefYZ += numNeighbors;
      rateOfDefZX += numNeighbors; rateOfDefZY += numNeighbors; rateOfDefZZ += numNeighbors;
      sigmaXXN += numNeighbors; sigmaXYN += numNeighbors; sigmaXZN += numNeighbors; 
      sigmaYXN += numNeighbors; sigmaYYN += numNeighbors; sigmaYZN += numNeighbors;
      sigmaZXN += numNeighbors; sigmaZYN += numNeighbors; sigmaZZN += numNeighbors;
      sigmaXXNP1 += numNeighbors; sigmaXYNP1 += numNeighbors; sigmaXZNP1 += numNeighbors; 
      sigmaYXNP1 += numNeighbors; sigmaYYNP1 += numNeighbors; sigmaYZNP1 += numNeighbors;
      sigmaZXNP1 += numNeighbors; sigmaZYNP1 += numNeighbors; sigmaZZNP1 += numNeighbors;
      vmStress += numNeighbors; 
    }
  }
}

// Explicit template instantiation for double
template void updateElasticCauchyStress<double>
(
    const double* deltaTemperatureN,
    const double* deltaTemperatureNP1,
    const double* unrotatedRateOfDeformation,
    const double* unrotatedCauchyStressN,
    double* unrotatedCauchyStressNP1,
    int numPoints,
    double bulkMod,
    double shearMod,
    double alpha,
    double dt
);


template void updateElasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* unrotatedCauchyStressN, 
    double* unrotatedCauchyStressNP1, 
    double* vonMisesStress,
    int numPoints, 
    double bulkMod,
    double shearMod,
    double dt
);

template void updateElasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* unrotatedCauchyStressN, 
    double* unrotatedCauchyStressNP1, 
    double* vonMisesStress,
    const double* flyingPointFlag,
    int numPoints, 
    double bulkMod,
    double shearMod,
    double dt
);

template void updateBondLevelElasticCauchyStress<double>
(
    const double* bondLevelUnrotatedRateOfDeformationXX, 
    const double* bondLevelUnrotatedRateOfDeformationXY, 
    const double* bondLevelUnrotatedRateOfDeformationXZ, 
    const double* bondLevelUnrotatedRateOfDeformationYX, 
    const double* bondLevelUnrotatedRateOfDeformationYY, 
    const double* bondLevelUnrotatedRateOfDeformationYZ, 
    const double* bondLevelUnrotatedRateOfDeformationZX, 
    const double* bondLevelUnrotatedRateOfDeformationZY, 
    const double* bondLevelUnrotatedRateOfDeformationZZ, 
    const double* bondLevelUnrotatedCauchyStressXXN, 
    const double* bondLevelUnrotatedCauchyStressXYN, 
    const double* bondLevelUnrotatedCauchyStressXZN, 
    const double* bondLevelUnrotatedCauchyStressYXN, 
    const double* bondLevelUnrotatedCauchyStressYYN, 
    const double* bondLevelUnrotatedCauchyStressYZN, 
    const double* bondLevelUnrotatedCauchyStressZXN, 
    const double* bondLevelUnrotatedCauchyStressZYN, 
    const double* bondLevelUnrotatedCauchyStressZZN, 
    double* bondLevelUnrotatedCauchyStressXXNP1, 
    double* bondLevelUnrotatedCauchyStressXYNP1, 
    double* bondLevelUnrotatedCauchyStressXZNP1, 
    double* bondLevelUnrotatedCauchyStressYXNP1, 
    double* bondLevelUnrotatedCauchyStressYYNP1, 
    double* bondLevelUnrotatedCauchyStressYZNP1, 
    double* bondLevelUnrotatedCauchyStressZXNP1, 
    double* bondLevelUnrotatedCauchyStressZYNP1, 
    double* bondLevelUnrotatedCauchyStressZZNP1, 
    double* bondLevelVonMisesStress,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
);

template void updateBondLevelElasticCauchyStress<double>
(
    const double* bondLevelUnrotatedRateOfDeformationXX, 
    const double* bondLevelUnrotatedRateOfDeformationXY, 
    const double* bondLevelUnrotatedRateOfDeformationXZ, 
    const double* bondLevelUnrotatedRateOfDeformationYX, 
    const double* bondLevelUnrotatedRateOfDeformationYY, 
    const double* bondLevelUnrotatedRateOfDeformationYZ, 
    const double* bondLevelUnrotatedRateOfDeformationZX, 
    const double* bondLevelUnrotatedRateOfDeformationZY, 
    const double* bondLevelUnrotatedRateOfDeformationZZ, 
    const double* bondLevelUnrotatedCauchyStressXXN, 
    const double* bondLevelUnrotatedCauchyStressXYN, 
    const double* bondLevelUnrotatedCauchyStressXZN, 
    const double* bondLevelUnrotatedCauchyStressYXN, 
    const double* bondLevelUnrotatedCauchyStressYYN, 
    const double* bondLevelUnrotatedCauchyStressYZN, 
    const double* bondLevelUnrotatedCauchyStressZXN, 
    const double* bondLevelUnrotatedCauchyStressZYN, 
    const double* bondLevelUnrotatedCauchyStressZZN, 
    double* bondLevelUnrotatedCauchyStressXXNP1, 
    double* bondLevelUnrotatedCauchyStressXYNP1, 
    double* bondLevelUnrotatedCauchyStressXZNP1, 
    double* bondLevelUnrotatedCauchyStressYXNP1, 
    double* bondLevelUnrotatedCauchyStressYYNP1, 
    double* bondLevelUnrotatedCauchyStressYZNP1, 
    double* bondLevelUnrotatedCauchyStressZXNP1, 
    double* bondLevelUnrotatedCauchyStressZYNP1, 
    double* bondLevelUnrotatedCauchyStressZZNP1, 
    double* bondLevelVonMisesStress,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */




template<typename ScalarT>
void updateElasticCauchyStressSmallDef
(
ScalarT* DeformationGradient, 
ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
int numPoints, 
const ScalarT Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental
)
{
  // Hooke's law
  const ScalarT* defGrad = DeformationGradient;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  
  ScalarT strain[9];
  ScalarT C[6][6];
  ScalarT rotationMat[3][3], rotationMatX[3][3], rotationMatY[3][3], rotationMatZ[3][3], temp[3][3];
  double alpha[3];

  // 0 -> xx,  1 -> xy, 2 -> xz
  // 3 -> yx,  4 -> yy, 5 -> yz
  // 6 -> zx,  7 -> zy, 8 -> zz


  for(int iID=0 ; iID<numPoints ; ++iID, 
        defGrad+=9, sigmaNP1+=9, sigmaN+=9){
            // df muss genutzt werden EQ 34 Silling Stability --> d.h. ich muss das Delta bestimmen
            // Schadensmodell muss angepasst werden
            // Referenzkonfiguration muss angepasst werden
            // volle C Matrix muss hier rein
      alpha[0] = angles[3*iID  ];
      alpha[1] = angles[3*iID+1];
      alpha[2] = angles[3*iID+2];
      //CORRESPONDENCE::createRotatedPythonBasedStiff(Cstiff,alpha,C); 
      //if (type != 0){
      //  // rotation mathematical positive sense (downwards from x-axis)
        CORRESPONDENCE::createRotationMatrix(angles[3*iID  ],rotationMatX,0);
        CORRESPONDENCE::createRotationMatrix(angles[3*iID+1],rotationMatY,1);
        CORRESPONDENCE::createRotationMatrix(angles[3*iID+2],rotationMatZ,2);
        CORRESPONDENCE::MatrixMultiply3x3(rotationMatX, rotationMatY, temp);
        CORRESPONDENCE::MatrixMultiply3x3(temp, rotationMatZ, rotationMat);
        CORRESPONDENCE::createRotatedStiff(Cstiff,rotationMat,C);
      //
      //}
      //else{
      //  C[0][0] = Cstiff[0][0];C[1][0] = Cstiff[1][0];C[2][0] = Cstiff[2][0];C[3][0] = Cstiff[3][0];C[4][0] = Cstiff[4][0];C[5][0] = Cstiff[5][0];
      //  C[0][1] = Cstiff[0][1];C[1][1] = Cstiff[1][1];C[2][1] = Cstiff[2][1];C[3][1] = Cstiff[3][1];C[4][1] = Cstiff[4][1];C[5][1] = Cstiff[5][1];
      //  C[0][2] = Cstiff[0][2];C[1][2] = Cstiff[1][2];C[2][2] = Cstiff[2][2];C[3][2] = Cstiff[3][2];C[4][2] = Cstiff[4][2];C[5][2] = Cstiff[5][2];
      //  C[0][3] = Cstiff[0][3];C[1][3] = Cstiff[1][3];C[2][3] = Cstiff[2][3];C[3][3] = Cstiff[3][3];C[4][3] = Cstiff[4][3];C[5][3] = Cstiff[5][3];
      //  C[0][4] = Cstiff[0][4];C[1][4] = Cstiff[1][4];C[2][4] = Cstiff[2][4];C[3][4] = Cstiff[3][4];C[4][4] = Cstiff[4][4];C[5][4] = Cstiff[5][4];
      //  C[0][5] = Cstiff[0][5];C[1][5] = Cstiff[1][5];C[2][5] = Cstiff[2][5];C[3][5] = Cstiff[3][5];C[4][5] = Cstiff[4][5];C[5][5] = Cstiff[5][5];
      //}
    if (type==0){
        strain[0] = 0.5 * ( *(defGrad)*   *(defGrad)   + *(defGrad+3)* *(defGrad+3) + *(defGrad+6) * *(defGrad+6)  - 1.0 );
        strain[1] = 0.5 * ( *(defGrad)*   *(defGrad+1) + *(defGrad+3)* *(defGrad+4) + *(defGrad+6) * *(defGrad+7)  );
        strain[2] = 0.5 * ( *(defGrad)*   *(defGrad+2) + *(defGrad+3)* *(defGrad+5) + *(defGrad+6) * *(defGrad+8)  );
        strain[3] = 0.5 * ( *(defGrad)*   *(defGrad+1) + *(defGrad+3)* *(defGrad+4) + *(defGrad+6) * *(defGrad+7)  );
        strain[4] = 0.5 * ( *(defGrad+1)* *(defGrad+1) + *(defGrad+4)* *(defGrad+4) + *(defGrad+7) * *(defGrad+7)  - 1.0 );
        strain[5] = 0.5 * ( *(defGrad+1)* *(defGrad+2) + *(defGrad+4)* *(defGrad+5) + *(defGrad+7) * *(defGrad+8)  );
        strain[6] = 0.5 * ( *(defGrad)*   *(defGrad+2) + *(defGrad+3)* *(defGrad+5) + *(defGrad+6) * *(defGrad+8)  );
        strain[7] = 0.5 * ( *(defGrad+1)* *(defGrad+2) + *(defGrad+4)* *(defGrad+5) + *(defGrad+7) * *(defGrad+8)  );
        strain[8] = 0.5 * ( *(defGrad+2)* *(defGrad+2) + *(defGrad+5)* *(defGrad+5) + *(defGrad+8) * *(defGrad+8)  - 1.0 );
        
        *(sigmaNP1)   = C[0][0]*strain[0] + C[0][1]*strain[4] + C[0][2]*strain[8] + C[0][3]*(strain[5]+strain[7]) + C[0][4]*(strain[2]+strain[6]) + C[0][5]*(strain[1]+strain[3]);
        *(sigmaNP1+1) = C[5][0]*strain[0] + C[5][1]*strain[4] + C[5][2]*strain[8] + C[5][3]*(strain[5]+strain[7]) + C[5][4]*(strain[2]+strain[6]) + C[5][5]*(strain[1]+strain[3]);
        *(sigmaNP1+2) = C[4][0]*strain[0] + C[4][1]*strain[4] + C[4][2]*strain[8] + C[4][3]*(strain[5]+strain[7]) + C[4][4]*(strain[2]+strain[6]) + C[4][5]*(strain[1]+strain[3]);
        *(sigmaNP1+3) = *(sigmaNP1+1);
        *(sigmaNP1+4) = C[1][0]*strain[0] + C[1][1]*strain[4] + C[1][2]*strain[8] + C[1][3]*(strain[5]+strain[7]) + C[1][4]*(strain[2]+strain[6]) + C[1][5]*(strain[1]+strain[3]);
        *(sigmaNP1+5) = C[3][0]*strain[0] + C[3][1]*strain[4] + C[3][2]*strain[8] + C[3][3]*(strain[5]+strain[7]) + C[3][4]*(strain[2]+strain[6]) + C[3][5]*(strain[1]+strain[3]);
        *(sigmaNP1+6) = *(sigmaNP1+2);
        *(sigmaNP1+7) = *(sigmaNP1+5);
        *(sigmaNP1+8) =  C[2][0]*strain[0] + C[2][1]*strain[4] + C[2][2]*strain[8] + C[2][3]*(strain[5]+strain[7]) + C[2][4]*(strain[2]+strain[6]) + C[2][5]*(strain[1]+strain[3]);
        
      }
      if (type!=0){

          strain[0] = 0.5 * ( *(defGrad)*   *(defGrad)   + *(defGrad+3)* *(defGrad+3) - 1.0 );
          strain[1] = 0.5 * ( *(defGrad)*   *(defGrad+1) + *(defGrad+3)* *(defGrad+4) );
          strain[3] = 0.5 * ( *(defGrad+1)* *(defGrad)   + *(defGrad+4)* *(defGrad+3) );
          strain[4] = 0.5 * ( *(defGrad+1)* *(defGrad+1) + *(defGrad+4)* *(defGrad+4) - 1.0 );
        *(sigmaNP1)   = C[0][0]*strain[0] + C[0][1]*strain[4] + C[0][5]*(strain[1]+strain[3]);
        *(sigmaNP1+1) = C[5][0]*strain[0] + C[5][1]*strain[4] + C[5][5]*(strain[1]+strain[3]);

        *(sigmaNP1+2) = 0.0;

        *(sigmaNP1+3) = *(sigmaNP1+1);
        *(sigmaNP1+4) = C[1][0] *strain[0] + C[1][1]* strain[4] + C[1][5]*(strain[1]+strain[3]);
        *(sigmaNP1+5) = 0.0;
        *(sigmaNP1+6) = 0.0;
        *(sigmaNP1+7) = 0.0;
        *(sigmaNP1+8) = 0.0;
        
      }
      //if (incremental == true){
      //   for (int i = 0; i < 9; i++) {
      //        *(sigmaNP1+i) += *(sigmaN+i);
      //    }
      //}
  }     
}       


// Explicit template instantiation for double
template void updateElasticCauchyStressSmallDef<double>
(
double* DeformationGradient, 
double* unrotatedCauchyStressN, 
double* unrotatedCauchyStressNP1, 
int numPoints, 
const double Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental
);
template void updateElasticCauchyStressSmallDef<Sacado::Fad::DFad<double> >
(
Sacado::Fad::DFad<double>* DeformationGradient, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressN, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1, 
int numPoints, 
const Sacado::Fad::DFad<double> Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental
);



/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template<typename ScalarT>
void createRotationMatrix
(
const double alpha,
ScalarT rotMat[][3],
const int check
){
    const double PI  = 3.141592653589793238463;
    double rad = alpha*(PI)/180;
    
    
    
    if (check == 0){// x - direction
        rotMat[0][0] = 1; rotMat[0][1] = 0;        rotMat[0][2] = 0;
        rotMat[1][0] = 0; rotMat[1][1] = cos(rad); rotMat[1][2] = -sin(rad);
        rotMat[2][0] = 0; rotMat[2][1] = sin(rad); rotMat[2][2] =  cos(rad);
    }
    if (check == 1){// y - direction
        rotMat[0][0] =  cos(rad); rotMat[0][1] = 0; rotMat[0][2] = sin(rad);
        rotMat[1][0] = 0;         rotMat[1][1] = 1; rotMat[1][2] = 0;
        rotMat[2][0] = -sin(rad); rotMat[2][1] = 0; rotMat[2][2] = cos(rad);
    }
    if (check == 2){// z - direction
        rotMat[0][0] = cos(rad); rotMat[0][1] = -sin(rad); rotMat[0][2] = 0;
        rotMat[1][0] = sin(rad); rotMat[1][1] =  cos(rad); rotMat[1][2] = 0;
        rotMat[2][0] = 0;        rotMat[2][1] = 0;         rotMat[2][2] = 1;

    }
}


//// Explicit template instantiation for double
template void createRotationMatrix<double>
(
const double alpha,
double rotMat[][3],
const int check
);
template void createRotationMatrix<Sacado::Fad::DFad<double> >
(
const double alpha,
Sacado::Fad::DFad<double> rotMat[][3],
const int check
);



/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template<typename ScalarT>
void createRotatedStiff
(
const ScalarT C[][6],
const ScalarT rotMat[][3],
ScalarT Cnew[][6]
){
    
     //  rM[0,0]**2     , rM[0,1]**2     , rM[0,2]**2     , 2.*rM[0,1]*rM[0,2]             , 2.*rM[0,0]*rM[0,2]             , 2.*rM[0,0]*rM[0,1]              ],
     //[ rM[1,0]**2     , rM[1,1]**2     , rM[1,2]**2     , 2.*rM[1,1]*rM[1,2]             , 2.*rM[1,0]*rM[1,2]             , 2.*rM[1,0]*rM[1,1]              ],
     //[ rM[2,0]**2     , rM[2,1]**2     , rM[2,2]**2     , 2.*rM[2,1]*rM[2,2]             , 2.*rM[2,0]*rM[2,2]             , 2.*rM[2,0]*rM[2,1]              ],
     //[ rM[1,0]*rM[2,0], rM[1,1]*rM[2,1], rM[1,2]*rM[2,2], rM[1,1]*rM[2,2]+rM[1,2]*rM[2,1], rM[1,0]*rM[2,2]+rM[1,2]*rM[2,0], rM[1,0]*rM[2,1]+rM[1,1]*rM[2,0] ],
     //[ rM[0,0]*rM[2,0], rM[0,1]*rM[2,1], rM[0,2]*rM[2,2], rM[0,1]*rM[2,2]+rM[0,2]*rM[2,1], rM[0,0]*rM[2,2]+rM[0,2]*rM[2,0], rM[0,0]*rM[2,1]+rM[0,1]*rM[2,0] ],
     //[ rM[0,0]*rM[1,0], rM[0,1]*rM[1,1], rM[0,2]*rM[1,2], rM[0,1]*rM[1,2]+rM[0,2]*rM[1,1], rM[0,0]*rM[1,2]+rM[0,2]*rM[1,0], rM[0,0]*rM[1,1]+rM[0,1]*rM[1,0] ],
     //])
    
    
    ScalarT tm[6][6];
    bool GlobToLoc = false; 
    if (GlobToLoc){
    
        tm[0][0] = rotMat[0][0]*rotMat[0][0]; tm[0][1] = rotMat[0][1]*rotMat[0][1]; tm[0][2] = rotMat[0][2]*rotMat[0][2]; tm[0][3] = 2*rotMat[0][1]*rotMat[0][2]; tm[0][4] = 2*rotMat[0][0]*rotMat[0][2]; tm[0][5] = 2*rotMat[0][0]*rotMat[0][1];
        tm[1][0] = rotMat[1][0]*rotMat[1][0]; tm[1][1] = rotMat[1][1]*rotMat[1][1]; tm[1][2] = rotMat[1][2]*rotMat[1][2]; tm[1][3] = 2*rotMat[1][1]*rotMat[1][2]; tm[1][4] = 2*rotMat[1][0]*rotMat[1][2]; tm[1][5] = 2*rotMat[1][0]*rotMat[1][1];  
        tm[2][0] = rotMat[2][0]*rotMat[2][0]; tm[2][1] = rotMat[2][1]*rotMat[2][1]; tm[2][2] = rotMat[2][2]*rotMat[2][2]; tm[2][3] = 2*rotMat[2][1]*rotMat[2][2]; tm[2][4] = 2*rotMat[2][0]*rotMat[2][2]; tm[2][5] = 2*rotMat[2][0]*rotMat[2][1];
        tm[3][0] = rotMat[1][0]*rotMat[2][0]; tm[3][1] = rotMat[1][1]*rotMat[2][1]; tm[3][2] = rotMat[1][2]*rotMat[2][2]; tm[3][3] = rotMat[1][1]*rotMat[2][2]+rotMat[1][2]*rotMat[2][1]; tm[3][4] = rotMat[1][0]*rotMat[2][2]+rotMat[1][2]*rotMat[2][0]; tm[3][5] =rotMat[1][0]*rotMat[2][1]+rotMat[1][1]*rotMat[2][0];
        tm[4][0] = rotMat[0][0]*rotMat[2][0]; tm[4][1] = rotMat[0][1]*rotMat[2][1]; tm[4][2] = rotMat[0][2]*rotMat[2][2]; tm[4][3] = rotMat[0][1]*rotMat[2][2]+rotMat[0][2]*rotMat[2][1]; tm[4][4] = rotMat[0][0]*rotMat[2][2]+rotMat[0][2]*rotMat[2][0]; tm[4][5] =rotMat[0][0]*rotMat[2][1]+rotMat[0][1]*rotMat[2][0];
        tm[5][0] = rotMat[0][0]*rotMat[1][0]; tm[5][1] = rotMat[0][1]*rotMat[1][1]; tm[5][2] = rotMat[0][2]*rotMat[1][2]; tm[5][3] = rotMat[0][1]*rotMat[1][2]+rotMat[0][2]*rotMat[1][1]; tm[5][4] = rotMat[0][0]*rotMat[1][2]+rotMat[0][2]*rotMat[1][0]; tm[5][5] =rotMat[0][0]*rotMat[1][1]+rotMat[0][1]*rotMat[1][0];
    }
    else{
        tm[0][0] = rotMat[0][0]*rotMat[0][0];   tm[0][1] = rotMat[0][1]*rotMat[0][1];   tm[0][2] = rotMat[0][2]*rotMat[0][2];   tm[0][3] = rotMat[0][1]*rotMat[0][2]; tm[0][4] = rotMat[0][0]*rotMat[0][2]; tm[0][5] = rotMat[0][0]*rotMat[0][1];
        tm[1][0] = rotMat[1][0]*rotMat[1][0];   tm[1][1] = rotMat[1][1]*rotMat[1][1];   tm[1][2] = rotMat[1][2]*rotMat[1][2];   tm[1][3] = rotMat[1][1]*rotMat[1][2]; tm[1][4] = rotMat[1][0]*rotMat[1][2]; tm[1][5] = rotMat[1][0]*rotMat[1][1];  
        tm[2][0] = rotMat[2][0]*rotMat[2][0];   tm[2][1] = rotMat[2][1]*rotMat[2][1];   tm[2][2] = rotMat[2][2]*rotMat[2][2];   tm[2][3] = rotMat[2][1]*rotMat[2][2]; tm[2][4] = rotMat[2][0]*rotMat[2][2]; tm[2][5] = rotMat[2][0]*rotMat[2][1];
        tm[3][0] = 2*rotMat[1][0]*rotMat[2][0]; tm[3][1] = 2*rotMat[1][1]*rotMat[2][1]; tm[3][2] = 2*rotMat[1][2]*rotMat[2][2]; tm[3][3] = rotMat[1][1]*rotMat[2][2]+rotMat[1][2]*rotMat[2][1]; tm[3][4] = rotMat[1][0]*rotMat[2][2]+rotMat[1][2]*rotMat[2][0]; tm[3][5] =rotMat[1][0]*rotMat[2][1]+rotMat[1][1]*rotMat[2][0];
        tm[4][0] = 2*rotMat[0][0]*rotMat[2][0]; tm[4][1] = 2*rotMat[0][1]*rotMat[2][1]; tm[4][2] = 2*rotMat[0][2]*rotMat[2][2]; tm[4][3] = rotMat[0][1]*rotMat[2][2]+rotMat[0][2]*rotMat[2][1]; tm[4][4] = rotMat[0][0]*rotMat[2][2]+rotMat[0][2]*rotMat[2][0]; tm[4][5] =rotMat[0][0]*rotMat[2][1]+rotMat[0][1]*rotMat[2][0];
        tm[5][0] = 2*rotMat[0][0]*rotMat[1][0]; tm[5][1] = 2*rotMat[0][1]*rotMat[1][1]; tm[5][2] = 2*rotMat[0][2]*rotMat[1][2]; tm[5][3] = rotMat[0][1]*rotMat[1][2]+rotMat[0][2]*rotMat[1][1]; tm[5][4] = rotMat[0][0]*rotMat[1][2]+rotMat[0][2]*rotMat[1][0]; tm[5][5] =rotMat[0][0]*rotMat[1][1]+rotMat[0][1]*rotMat[1][0];
        
        
    }
    ScalarT Ctemp[6][6];
    bool transpose = true;
    CORRESPONDENCE::MatMul(6,tm ,C ,Ctemp, transpose);
    transpose = false;
    CORRESPONDENCE::MatMul(6,Ctemp ,tm ,Cnew, transpose);
}
template void createRotatedStiff<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double> C[][6],
const Sacado::Fad::DFad<double> rotMat[][3],
Sacado::Fad::DFad<double>  Cnew[][6]
);


// Explicit template instantiation for double
template void createRotatedStiff<double>
(
const double C[][6],
const double rotationMat[][3],
double Cnew[][6]
);






}


