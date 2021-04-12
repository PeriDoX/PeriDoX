/*! \file Peridigm_EnergyReleaseDamageCorrepondenceModel.cpp */

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
// Author of this Routine
// Christian Willberg   christian.willberg@dlr.de
// German Aerospace Center
//@HEADER

#include "Peridigm_EnergyReleaseDamageCorrepondenceModel.hpp"
#include "Peridigm_Field.hpp"
#include "material_utilities.h"
#include "correspondence.h"
#include <thread>
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::EnergyReleaseDamageCorrepondenceModel::EnergyReleaseDamageCorrepondenceModel(const Teuchos::ParameterList& params)
: DamageModel(params), 
m_applyThermalStrains(false),
m_modelCoordinatesFieldId(-1),
m_coordinatesFieldId(-1),
m_damageFieldId(-1),
m_bondDamageFieldId(-1),
m_deltaTemperatureFieldId(-1),
m_dilatationFieldId(-1),
m_weightedVolumeFieldId(-1),
m_horizonFieldId(-1),
m_piolaStressTimesInvShapeTensorXId(-1),
m_piolaStressTimesInvShapeTensorYId(-1),
m_piolaStressTimesInvShapeTensorZId(-1),
m_detachedNodesFieldId(-1),
m_forceDensityFieldId(-1),
m_deformationGradientFieldId(-1),
m_hourglassStiffId(-1),
m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()) {

    
    if (params.isParameter("Critical Energy")) {
        m_criticalEnergyTension = params.get<double>("Critical Energy");
          
    } 
    if (params.isParameter("Degradation Factor")){
        degradationFactor = 1;
        degradationFactor = params.get<double>("Degradation Factor");
        
    }else{
        degradationFactor = 1; 
    }
    m_incremental = true;
    if (params.isParameter("Incremental")){
       m_incremental = params.get<bool>("Incremental");
    }
    m_criticalEnergyTension = params.get<double>("Critical Energy");
    m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    m_planeStrain=false;
    m_planeStress=false;
    detachedNodesCheck = false;
    if (params.isParameter("Detached Nodes Check"))
        detachedNodesCheck = params.get<bool>("Detached Nodes Check");
    if (params.isParameter("Plane Strain"))
      m_planeStrain = params.get<bool>("Plane Strain");
    if (params.isParameter("Plane Stress"))
      m_planeStress = params.get<bool>("Plane Stress");
    m_plane = false;
   
    m_onlyTension = false;
    if(params.isParameter("Only Tension")){
        m_onlyTension = params.get<bool>("Only Tension");

    }

  //************************************
  // wie komme ich an den Namen??
  //************************************
  //if (params.isParameter("Linear Elastic Correspondence")){
   
    //std::cout<<"Use Material: Linear Elastic Correspondence"<<std::endl;

    if (m_planeStrain==true){
        //m_plane=true;
        m_plane = true;
        m_Thickness = params.get<double>("Thickness");
        //std::cout<<"Method 2D Plane Strain"<<std::endl;
    }
    if (m_planeStress==true){
       // m_plane=true;
        m_plane = true;
        m_Thickness = params.get<double>("Thickness");
        //std::cout<<"WRN: Method 2D Plane Stress --> not fully implemented yet"<<std::endl;
    }
   
// Params anpassen. Irgendwo muss das befuellt werden, da bspw. die thermische Verformung bei Material
// und Schaden existiert
//
//

    m_pi = 3.14159;

    if (params.isParameter("Thermal Expansion Coefficient")) {
        m_alpha = params.get<double>("Thermal Expansion Coefficient");
        m_applyThermalStrains = true;
    }

    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT,"Model_Coordinates");
    m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
    m_dilatationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
    m_damageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
    m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
    m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
    m_piolaStressTimesInvShapeTensorXId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorX");
    m_piolaStressTimesInvShapeTensorYId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorY");
    m_piolaStressTimesInvShapeTensorZId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorZ");
    m_detachedNodesFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
    m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
    m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
    m_hourglassStiffId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Hourglass_Stiffness");
 
    
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_modelCoordinatesFieldId);
    m_fieldIds.push_back(m_coordinatesFieldId);
    m_fieldIds.push_back(m_weightedVolumeFieldId);
    m_fieldIds.push_back(m_dilatationFieldId);
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_damageFieldId);
    m_fieldIds.push_back(m_detachedNodesFieldId);
    m_fieldIds.push_back(m_bondDamageFieldId);
    m_fieldIds.push_back(m_horizonFieldId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorXId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorYId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorZId);
    m_fieldIds.push_back(m_forceDensityFieldId);
    m_fieldIds.push_back(m_hourglassStiffId);
    m_fieldIds.push_back(m_deformationGradientFieldId);

}

PeridigmNS::EnergyReleaseDamageCorrepondenceModel::~EnergyReleaseDamageCorrepondenceModel() {
}

void
PeridigmNS::EnergyReleaseDamageCorrepondenceModel::initialize(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager) const {
    double *damage, *bondDamage;


    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
    
    dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    // Initialize damage to zero
    int neighborhoodListIndex(0);
    int bondIndex(0);
    int nodeId, numNeighbors;
    int iID, iNID;

    for (iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        damage[nodeId] = 0.0;
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        neighborhoodListIndex += numNeighbors;
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            bondDamage[bondIndex] = 0.0;
            bondIndex += 1;
        }
    }

}

void
PeridigmNS::EnergyReleaseDamageCorrepondenceModel::computeDamage(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager) const {

    double *x, *y, *damage, *bondDamageNP1, *horizon, *vol, *detachedNodes;
    
    
    double criticalEnergyTension(-1.0);
    // for temperature dependencies easy to extent
    double *deltaTemperature = NULL;
    double *tempStressX, *tempStressY, *tempStressZ;
    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);

    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);

    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&vol);

    dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);


    
    dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->ExtractView(&tempStressX);
    dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->ExtractView(&tempStressY);
    dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->ExtractView(&tempStressZ);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
    ////////////////////////////////////////////////////////
    // Hourglass control
    // not synchronized yet; Hourglass correction is used only from nodeId site
    //
    double *defGrad, *hourglassStiff;
    vector<double> TSvector(3), hourglassStiffVector(9);
    double* hStiff = &hourglassStiffVector[0];
    double* TS = &TSvector[0];
    ///////////////////////////////////////////////////////
    dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->ExtractView(&hourglassStiff);
    dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&defGrad);
    //std::cout<< "heredam"<<std::endl;
    // Set the bond damage to the previous value --> needed for iteration in implicit time integration
    *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));
     double *forceDensity;
    dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);
    
    ////////////////////////////////////////////////////
    double trialDamage(0.0);
    int neighborhoodListIndex(0), bondIndex(0);
    int nodeId, numNeighbors, neighborID, iID, iNID;
    double totalDamage;
    double nodeInitialX[3], nodeCurrentX[3];
    double omegaP1, omegaP2;
    double critIso;
    
    double quadhorizon;
    // bond force density state
    double TX, TY, TZ, TXN, TYN, TZN;
    // bond force density state projected to deformed bond
    double TPX, TPY, TPZ, TPXN, TPYN, TPZN;
    double factor, factorN;
    // deformed state
    double Y_dx, Y_dy, Y_dz;
    // initial state
    double X_dx, X_dy, X_dz;
    // displacement vector state and abolute value
    double etaX, etaY, etaZ, dEta, normEtaSq;
    double bondEnergy;
    double dX, dY;//, dYSq;

    //---------------------------
    // INITIALIZE PROCESS STEP t
    //---------------------------
    // Foster 2009 Journal for Multiscale Computational Engineering
    // m_criticalEnergy = 4 G / (pi delta^4)
    // what if two horizons meet? two criterions will hit the same bond..
    //double quadhorizon =  4 /( m_pi * m_horizon * m_horizon * m_horizon * m_horizon );
    
    if (m_criticalEnergyTension > 0.0)
        criticalEnergyTension = m_criticalEnergyTension;
 
    // Update the bond damage
    // Break bonds if the bond energy potential is greater than the critical bond energy potential
    //---------------------------
    // DAMAGE ANALYSIS
    //---------------------------
    bondIndex = 0;

    for (iID = 0; iID < numOwnedPoints; ++iID, defGrad+=9, hourglassStiff+=9) {
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        
        nodeId = ownedIDs[iID];
        nodeInitialX[0] = x[nodeId*3];
        nodeInitialX[1] = x[nodeId*3+1];
        nodeInitialX[2] = x[nodeId*3+2];
        nodeCurrentX[0] = y[nodeId*3];
        nodeCurrentX[1] = y[nodeId*3+1];
        nodeCurrentX[2] = y[nodeId*3+2];
        
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            
            
            neighborID = neighborhoodList[neighborhoodListIndex++];
            if (detachedNodes[nodeId]!=0) continue;
            if (detachedNodes[neighborID]!=0) continue;
            X_dx = x[neighborID*3]   - nodeInitialX[0];
            X_dy = x[neighborID*3+1] - nodeInitialX[1];
            X_dz = x[neighborID*3+2] - nodeInitialX[2];
 
            Y_dx = y[neighborID*3]   - nodeCurrentX[0];
            Y_dy = y[neighborID*3+1] - nodeCurrentX[1];
            Y_dz = y[neighborID*3+2] - nodeCurrentX[2];
            etaX  = (Y_dx-X_dx);
            etaY  = (Y_dy-X_dy);
            etaZ  = (Y_dz-X_dz);
            
                         
            //double uNx = y[neighborID*3]   - x[neighborID*3];
            //double uNy = y[neighborID*3+1] - x[neighborID*3+1];
            //double uNz = y[neighborID*3+2] - x[neighborID*3+2];
            
            
            dX = distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2], x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
            dY = distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2], y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
            dEta = dY-dX;
            

            normEtaSq = etaX*etaX+etaY*etaY+etaZ*etaZ;

            bool Tension = true;
            
            if (m_onlyTension == true && dEta<0) Tension = false;
            if (Tension == true){
                if (normEtaSq>0){
                    omegaP1 = MATERIAL_EVALUATION::scalarInfluenceFunction(dX, horizon[nodeId]); 
                    omegaP2 = MATERIAL_EVALUATION::scalarInfluenceFunction(-dX, horizon[neighborID]); 
                    // average Force has to be taken
                    // if not the case where 
                    if (m_plane == true) X_dz = 0;
                    double temp = 0.5*(1-bondDamageNP1[bondIndex]);
                    
                    double FxsiX = *(defGrad)   * X_dx + *(defGrad+1) * X_dy + *(defGrad+2) * X_dz;
                    double FxsiY = *(defGrad+3) * X_dx + *(defGrad+4) * X_dy + *(defGrad+5) * X_dz;
                    double FxsiZ = *(defGrad+6) * X_dx + *(defGrad+7) * X_dy + *(defGrad+8) * X_dz;
                    hStiff[0] = *(hourglassStiff  ); hStiff[1] = *(hourglassStiff+1);  hStiff[2] = *(hourglassStiff+2);
                    hStiff[3] = *(hourglassStiff+3); hStiff[4] = *(hourglassStiff+4);  hStiff[5] = *(hourglassStiff+5);
                    hStiff[6] = *(hourglassStiff+6); hStiff[7] = *(hourglassStiff+7);  hStiff[8] = *(hourglassStiff+8);
                    CORRESPONDENCE::computeCorrespondenceStabilityWanEtAlShort(FxsiX,FxsiY,FxsiZ,Y_dx,Y_dy,Y_dz,hStiff,TS);
               // std::cout<< TS[0]<<" dam"<< TS[1]<<std::endl;
                    // volumes, or volume relation to include? should be inside via the shape tensor, which is included in tempStress
                    TX  =   omegaP1 * ( tempStressX[3*nodeId]     * X_dx + tempStressX[3*nodeId+1]     * X_dy + tempStressX[3*nodeId+2]     * X_dz + m_hourglassCoefficient*TS[0]);
                    TY  =   omegaP1 * ( tempStressY[3*nodeId]     * X_dx + tempStressY[3*nodeId+1]     * X_dy + tempStressY[3*nodeId+2]     * X_dz + m_hourglassCoefficient*TS[1]);
                    TZ  =   omegaP1 * ( tempStressZ[3*nodeId]     * X_dx + tempStressZ[3*nodeId+1]     * X_dy + tempStressZ[3*nodeId+2]     * X_dz + m_hourglassCoefficient*TS[2]);
                    // undeformedBondX, undeformedBondY, undeformedBondZ of bond 1-2 equal to -undeformedBondX, -undeformedBondY, -undeformedBondZ of bond 2-1
                    TXN =   omegaP2 * ( tempStressX[3*neighborID] * X_dx + tempStressX[3*neighborID+1] * X_dy + tempStressX[3*neighborID+2] * X_dz + m_hourglassCoefficient*TS[0]);
                    TYN =   omegaP2 * ( tempStressY[3*neighborID] * X_dx + tempStressY[3*neighborID+1] * X_dy + tempStressY[3*neighborID+2] * X_dz + m_hourglassCoefficient*TS[1]);
                    TZN =   omegaP2 * ( tempStressZ[3*neighborID] * X_dx + tempStressZ[3*neighborID+1] * X_dy + tempStressZ[3*neighborID+2] * X_dz + m_hourglassCoefficient*TS[2]);
                    //std::cout<< "here2"<<std::endl;
                    // orthogonal projection of T and TN to the relative displacement vector Foster et al. "An energy based .."
                    // --> die senkrecht zur Projektion stehenden Anteile entsprechen eventuell den Schubanteilen. D.h. man könnte das Kriterium hier splitten.
                    factor = (etaX*TX + etaY*TY + etaZ*TZ)/normEtaSq;
                    TPX = factor*etaX; TPY = factor*etaY; TPZ = factor*etaZ;
                    
                    //factor = (Y_dx*TX + Y_dy*TY + Y_dz*TZ)/dYSq;
                    //TPX = factor*Y_dx; TPY = factor*Y_dy; TPZ = factor*Y_dz;
                    
                    factorN = (etaX*TXN + etaY*TYN + etaZ*TZN)/normEtaSq;
                    TPXN = factorN*etaX; TPYN = factorN*etaY; TPZN = factorN*etaZ;
                    
                    //factorN = (Y_dx*TXN + Y_dy*TYN + Y_dz*TZN)/dYSq;
                    //TPXN = factorN*Y_dx; TPYN = factorN*Y_dy; TPZN = factorN*Y_dz;
                    
                    bondEnergy = 0.5*temp*(abs(TPX*etaX)+abs(TPXN*etaX)+abs(TPY*etaY)+abs(TPYN*etaY)+abs(TPZ*etaZ)+abs(TPZN*etaZ)) ;
                    //if (iID == 1) std::cout<<bondEnergy<<std::endl;
                }
                else
                {
                    bondEnergy = 0;
                }
                //bondEnergy = 0.5*(sqrt(TX*TX)+sqrt(TXN*TXN))*sqrt((uNx-ux)*(uNx-ux)) + 0.5*(sqrt(TY*TY)+sqrt(TYN*TYN))*sqrt((uNy-uy)*(uNy-uy)) + 0.5*(sqrt(TZ*TZ)+sqrt(TZN*TZN))*sqrt((uNz-uz)*(uNz-uz));
                

                double avgHorizon = 0.5*(horizon[nodeId]+horizon[neighborID]);
                if (bondEnergy<0){
                    std::cout<<TPX<<" "<< TPXN<<" BE "<<bondEnergy<<std::endl;
                }
                ////////////////////////////////////////////////////////////////
                //--> to check, depth is not included yet. How to handle??
                ////////////////////////////////////////////////////////////////
                if (m_planeStrain==false&&m_planeStress==false){
                   quadhorizon =  4 /( m_pi * avgHorizon * avgHorizon * avgHorizon * avgHorizon );
                }
                else
                {

                   quadhorizon =  3 /( avgHorizon * avgHorizon * avgHorizon * m_Thickness );
                }
                //quadhorizon =  4 /( m_pi * avgHorizon * avgHorizon * avgHorizon * avgHorizon ); 
                critIso = bondEnergy/(criticalEnergyTension*quadhorizon);
                //critIso = 0;
                
                //critIso = bondEnergy/(criticalEnergyTension*quadhorizon);
                
               // std::cout<<bondEnergy<<" EE "<<criticalEnergyTension*quadhorizon<<" Stress1 "<<dEta<<std::endl;
                trialDamage = 0.0;
                if (criticalEnergyTension > 0.0 && critIso > 1.0) {
                    trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }

                if (trialDamage > bondDamageNP1[bondIndex]) {
                    if (trialDamage>1)trialDamage = 1;
                    bondDamageNP1[bondIndex] = trialDamage;

                }
            }
            bondIndex += 1;

            }

        }
    //  Update the element damage (percent of bonds broken)
    if (detachedNodesCheck == true){
        int check = 1;  //set check = 1 to start the loop
        //std::cout<< "detached"<<std::endl;
        while (check != 0){
            check = checkDetachedNodes(numOwnedPoints, ownedIDs, neighborhoodList, dataManager);
        }
    }
    neighborhoodListIndex = 0;
    bondIndex = 0;
    double volume;
    for (iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        //neighborhoodListIndex += numNeighbors;
        totalDamage = 0.0;
        volume = vol[nodeId];
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            
            neighborID = neighborhoodList[neighborhoodListIndex++];
            // must be zero to avoid synchronization errors
            
            totalDamage += bondDamageNP1[bondIndex]*vol[neighborID];
            volume += vol[neighborID];
            bondIndex += 1;
        }
        if (numNeighbors > 0)
            totalDamage /= numNeighbors;
        else
            totalDamage = 0.0;

        damage[nodeId] = totalDamage/volume;

    }
}

int PeridigmNS::EnergyReleaseDamageCorrepondenceModel::checkDetachedNodes(
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,                                                     
                                                      PeridigmNS::DataManager& dataManager
                                                      ) const
{
  double *bondDamage, *detachedNodes, *x, *volume, *horizon;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
 
  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;
  double nodeInitialX[3],  initialDistance;
  double checkShapeTensor[9], neighborVolume, checkShapeTensorInv[9];
  double omega;
  double determinant;
  int check = 0, matrixInversionReturnCode = 0;

  for(iID=0 ; iID<numOwnedPoints ; ++iID){
      nodeId = ownedIDs[iID];
      nodeInitialX[0] = x[nodeId*3];
      nodeInitialX[1] = x[nodeId*3+1];
      nodeInitialX[2] = x[nodeId*3+2];

      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      
      checkShapeTensor[0] = 0.0;
      checkShapeTensor[1] = 0.0;
      checkShapeTensor[2] = 0.0;
      checkShapeTensor[3] = 0.0;
      checkShapeTensor[4] = 0.0;
      checkShapeTensor[5] = 0.0;
      checkShapeTensor[6] = 0.0;
      checkShapeTensor[7] = 0.0;
      checkShapeTensor[8] = 0.0;
  
      for(iNID=0 ; iNID<numNeighbors ; ++iNID){
          
          neighborID = neighborhoodList[neighborhoodListIndex++];
          
          if (detachedNodes[nodeId]==0){
          
              neighborVolume = volume[neighborID];
              initialDistance = distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
                                  x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
              omega = MATERIAL_EVALUATION::scalarInfluenceFunction(initialDistance, horizon[iID]);
              //double omega = 1.0;
              double temp = (1.0 - bondDamage[bondIndex]) * omega * neighborVolume;
              
              double undeformedBondX =  x[neighborID*3]   - nodeInitialX[0];
              double undeformedBondY =  x[neighborID*3+1] - nodeInitialX[1];
              double undeformedBondZ =  x[neighborID*3+2] - nodeInitialX[2];
              
              checkShapeTensor[0]  += temp * undeformedBondX * undeformedBondX;
              checkShapeTensor[1]  += temp * undeformedBondX * undeformedBondY;
              checkShapeTensor[2]  += temp * undeformedBondX * undeformedBondZ;
              checkShapeTensor[3]  += temp * undeformedBondY * undeformedBondX;
              checkShapeTensor[4]  += temp * undeformedBondY * undeformedBondY;
              checkShapeTensor[5]  += temp * undeformedBondY * undeformedBondZ;
              checkShapeTensor[6]  += temp * undeformedBondZ * undeformedBondX;
              checkShapeTensor[7]  += temp * undeformedBondZ * undeformedBondY;
              checkShapeTensor[8]  += temp * undeformedBondZ * undeformedBondZ;
          }
          if (detachedNodes[neighborID]!=0&&bondDamage[bondIndex] != 1.0){// bondDamage check to avoid infinite loop; all bonds then destroyed
              bondDamage[bondIndex] = 1.0;
              check = 1;
          }
          
          bondIndex += 1;
      }
      
      if (detachedNodes[nodeId]!=0) continue;
      
      if (m_plane==true){
        matrixInversionReturnCode =
        CORRESPONDENCE::Invert2by2Matrix(checkShapeTensor, determinant, checkShapeTensorInv);
        }
      else{
        matrixInversionReturnCode =
        CORRESPONDENCE::Invert3by3Matrix(checkShapeTensor, determinant, checkShapeTensorInv);
        }
      
      if (matrixInversionReturnCode != 0){// to be checked
          check = 1;
          detachedNodes[nodeId]=1.;
          int bondIndex2 = bondIndex-numNeighbors; // set index back, to have the same correct bonds
          // delete all connected bonds
          for(iNID=0 ; iNID<numNeighbors ; ++iNID){
                  bondDamage[bondIndex2] = 1.0;
                  bondIndex2 += 1;
          }   
        matrixInversionReturnCode = 0;
      }
  }

  return check;
}


