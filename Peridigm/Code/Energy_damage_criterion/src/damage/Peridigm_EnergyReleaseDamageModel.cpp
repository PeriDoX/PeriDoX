/*! \file Peridigm_EnergyReleaseDamageModel.cpp */

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
//////////////////////////////////////////////////////////////////////////////////////    
// Routine developed for Peridigm by
// DLR Composite Structures and Adaptive Systems
//                                __/|__
//                                /_/_/_/  
//            www.dlr.de/fa/en      |/ DLR
//////////////////////////////////////////////////////////////////////////////////////    
// Questions?
// Christian Willberg  christian.willberg@dlr.de
// Martin Raedel       martin.raedel@dlr.de
//////////////////////////////////////////////////////////////////////////////////////   
//@HEADER

#include "Peridigm_EnergyReleaseDamageModel.hpp"
#include "Peridigm_Field.hpp"
#include "material_utilities.h"
#include <thread>
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::EnergyReleaseDamageModel::EnergyReleaseDamageModel(const Teuchos::ParameterList& params)
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
m_damageModelFieldId(-1),
m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()) {

   
    if (params.isParameter("Critical Energy")) {
        m_criticalEnergyTension = params.get<double>("Critical Energy Tension");
        m_criticalEnergyCompression = params.get<double>("Critical Energy Compression");
        m_criticalEnergyShear = params.get<double>("Critical Energy Shear");
          
    } else {
        if (params.isParameter("Critical Energy Tension"))
            m_criticalEnergyTension = params.get<double>("Critical Energy Tension");
        else
            m_criticalEnergyTension = -1.0;
        if (params.isParameter("Critical Energy Compression"))
            m_criticalEnergyCompression = params.get<double>("Critical Energy Compression");
        else
            m_criticalEnergyCompression = -1.0;
        if (params.isParameter("Critical Energy Shear"))
            m_criticalEnergyShear = params.get<double>("Critical Energy Shear");
        else
            m_criticalEnergyShear = -1.0; 
        m_type = 1;
        if (params.isType<string>("Energy Criterion"))
            m_type = 1;
        if (params.isType<string>("Power Law"))
            m_type = 2;
        if (params.isType<string>("Separated"))
            m_type = 3; 
    }

    m_pi = 3.14159;

    if (params.isParameter("Thermal Expansion Coefficient")) {
        m_alpha = params.get<double>("Thermal Expansion Coefficient");
        m_applyThermalStrains = true;
    }

    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
    m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
    m_volumeFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
    m_dilatationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
    m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
    m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
    m_horizonFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::CONSTANT, "Horizon");
    m_damageModelFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE, PeridigmNS::PeridigmField::VECTOR, PeridigmNS::PeridigmField::TWO_STEP, "Damage_Model_Data");
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_modelCoordinatesFieldId);
    m_fieldIds.push_back(m_coordinatesFieldId);
    m_fieldIds.push_back(m_weightedVolumeFieldId);
    m_fieldIds.push_back(m_dilatationFieldId);
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_damageFieldId);
    m_fieldIds.push_back(m_damageModelFieldId);
    m_fieldIds.push_back(m_bondDamageFieldId);
    m_fieldIds.push_back(m_horizonFieldId);

}

PeridigmNS::EnergyReleaseDamageModel::~EnergyReleaseDamageModel() {
}

void
PeridigmNS::EnergyReleaseDamageModel::initialize(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager) const {
    double *damage, *bondDamage;


    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);


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
PeridigmNS::EnergyReleaseDamageModel::computeDamage(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager) const {

    double *x, *y, *damage, *bondDamageNP1, *horizon;
    double *cellVolume, *weightedVolume, *damageModel;
    double criticalEnergyTension(-1.0), criticalEnergyCompression(-1.0), criticalEnergyShear(-1.0);
    // for temperature dependencies easy to extent
    double *deltaTemperature = NULL;
    double m_alpha = 0;

    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
    //////////////////////////////////////////////////////////////
    // transfer of data is done in ComputeDilation --> PeridigmMaterial.cpp
    //////////////////////////////////////////////////////////////
    dataManager.getData(m_damageModelFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageModel);
    dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    ////////////////////////////////////
    ////////////////////////////////////
    // transfer of data is done in ComputeDilation --> PeridigmMaterial.cpp
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
    

    ////////////////////////////////////////////////////
    double trialDamage(0.0);
    int neighborhoodListIndex(0), bondIndex(0);
    int nodeId, numNeighbors, neighborID, iID, iNID;
    double totalDamage;
    double alphaP1, alphaP2;
    double nodeInitialX[3], nodeCurrentX[3], relativeExtension(0.0);
    double bondEnergyIsotropic(0.0), bondEnergyDeviatoric(0.0);
    double omegaP1, omegaP2;
    double critDev, critIso, critComp;
    double BulkMod1, BulkMod2;
    double degradationFactor = 1.0; // Optional parameter if bond should be degradated and not fully destroyed instantaneously 
    double avgHorizon, quadhorizon;
    //---------------------------
    // INITIALIZE PROCESS STEP t
    //---------------------------
    
    if (m_criticalEnergyTension > 0.0)
        criticalEnergyTension = m_criticalEnergyTension;
    if (m_criticalEnergyCompression > 0.0)
        criticalEnergyCompression = m_criticalEnergyCompression;
    if (m_criticalEnergyShear > 0.0)
        criticalEnergyShear = m_criticalEnergyShear;

    // Update the bond damage
    // Break bonds if the bond energy potential is greater than the critical bond energy potential
    //---------------------------
    // DAMAGE ANALYSIS
    //---------------------------
    bondIndex = 0;
    
    for (iID = 0; iID < numOwnedPoints; ++iID) {
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        
        nodeId = ownedIDs[iID];
        nodeInitialX[0] = x[nodeId*3];
        nodeInitialX[1] = x[nodeId*3+1];
        nodeInitialX[2] = x[nodeId*3+2];
        nodeCurrentX[0] = y[nodeId*3];
        nodeCurrentX[1] = y[nodeId*3+1];
        nodeCurrentX[2] = y[nodeId*3+2];
 
        double dilatationP1 = damageModel[3*nodeId];
        
        alphaP1 = 15.0 * damageModel[3*nodeId+2]; // weightedVolume is already included in --> Perdigm_Material.cpp;
        
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            
            trialDamage = 0.0;
            neighborID = neighborhoodList[neighborhoodListIndex++];
            
            double zeta = 
            distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
                 x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);

            double dY = 
            distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
                 y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
 
            relativeExtension = (dY - zeta)/zeta;    
             // the direction switches between both bonds. This results in a switch of the forces
             // as well. Therefore, all forces and bond deformations are normalized.
            double eP = dY - zeta;
            //double normZeta = sqrt(zeta*zeta);
            
            alphaP2 = 15.0 * damageModel[3*neighborID+2]; // weightedVolume is already included in --> Perdigm_Material.cpp;
            
            BulkMod1 = damageModel[3*nodeId+1];     // weightedVolume is already included in the variable --> Perdigm_Material.cpp;
            BulkMod2 = damageModel[3*neighborID+1]; // weightedVolume is already included in the variable --> Perdigm_Material.cpp;
            
            
            double dilatationP2 = damageModel[3*neighborID];
            
            omegaP1 = MATERIAL_EVALUATION::scalarInfluenceFunction(zeta, horizon[nodeId]); 
            omegaP2 = MATERIAL_EVALUATION::scalarInfluenceFunction(zeta, horizon[neighborID]); 

            double eiP1 = dilatationP1 * zeta / 3.0;
            double eiP2 = dilatationP2 * zeta / 3.0;
            double tiP1 = 3*BulkMod1*omegaP1*dilatationP1*zeta;
            double tiP2 = 3*BulkMod2*omegaP2*dilatationP2*zeta;
            

            double edP1 = eP - sqrt(eiP1*eiP1);
            double edP2 = eP - sqrt(eiP2*eiP2);
            double tdP1 = omegaP1*alphaP1*edP1;
            double tdP2 = omegaP2*alphaP2*edP2;
            // absolute values of energy, because it is positive and coordinate errors are avoided.
            bondEnergyIsotropic  = (1.0 - bondDamageNP1[bondIndex])*(sqrt(tiP1*eiP1*tiP1*eiP1) + sqrt(tiP2*eiP2*tiP2*eiP2)); 
            bondEnergyDeviatoric = (1.0 - bondDamageNP1[bondIndex])*(sqrt(tdP1*edP1*tdP2*edP2) + sqrt(tdP2*edP2*tdP2*edP2));
            // the average horizon is taken, if multiple horizons are used
            avgHorizon = 0.5*(horizon[nodeId]+horizon[neighborID]);

            // geometrical part of the energy value
            quadhorizon =  16.0 /( m_pi * avgHorizon * avgHorizon * avgHorizon * avgHorizon );
            // the factor 16 can be split in three parts:
            // 4 comes from the integration from Foster et al. (2009) Journal for Multiscale Computational Engineering; 
            // 2 comes from the force split motivated by the bond based formulation Bobaru et al. (2017) "Handbook of Peridynamik Modeling", page 48 ; 
            // 2 comes from the energy formulation itself
            critIso = 0.0;
            if (relativeExtension>0&&criticalEnergyTension != -1.0){
               critIso = (bondEnergyIsotropic/(criticalEnergyTension*quadhorizon));
            }
            critDev = 0.0;
            if (criticalEnergyShear != -1.0){
               critDev = (bondEnergyDeviatoric/(criticalEnergyShear*quadhorizon));
            }
            critComp = 0.0;
            if (relativeExtension < 0.0 && criticalEnergyCompression != -1.0){
                critComp = (bondEnergyIsotropic/(criticalEnergyCompression*quadhorizon));
                critIso = 0.0;
            }
            if (m_type == 1){ // Energy Criterion by Foster et al.(2009) Journal for Multiscale Computational Engineering;
                double bondEnergy = sqrt((tdP1+tiP1)*(tdP1+tiP1)*eP*eP) + sqrt((tdP2 + tiP2)*(tdP2 + tiP2)*eP*eP);
                critIso = bondEnergy/(criticalEnergyTension*quadhorizon);
                if (m_criticalEnergyTension > 0.0 && critIso > 1.0 ) {
                    trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                } 
            }
            
            if (m_type == 2){// Power Law
    
                if (m_criticalEnergyTension > 0.0 &&critIso*critIso + critDev*critDev + critComp*critComp > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
            }
            if (m_type == 3){ // Separated
                if (m_criticalEnergyTension > 0.0  &&critIso > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
                if (m_criticalEnergyTension > 0.0 && critDev > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                } 
                if (m_criticalEnergyCompression > 0.0  && critComp > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
            }
      
            if (trialDamage > bondDamageNP1[bondIndex]) {
                if (trialDamage>1)trialDamage = 1;
                bondDamageNP1[bondIndex] = trialDamage;

            }

            bondIndex += 1;

            }

        }
    

    //  Update the element damage (percent of bonds broken)

    neighborhoodListIndex = 0;
    bondIndex = 0;
    for (iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        //neighborhoodListIndex += numNeighbors;
        damageModel[3*nodeId] = 0.0;
        damageModel[3*nodeId+1] = 0.0;
        damageModel[3*nodeId+2] = 0.0;
        totalDamage = 0.0;
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            
            neighborID = neighborhoodList[neighborhoodListIndex++];
            // must be zero to avoid synchronization errors
            damageModel[3*neighborID] = 0.0;
            damageModel[3*neighborID+1] = 0.0;
            damageModel[3*neighborID+2] = 0.0;
            
            totalDamage += bondDamageNP1[bondIndex];

            bondIndex += 1;
        }
        if (numNeighbors > 0)
            totalDamage /= numNeighbors;
        else
            totalDamage = 0.0;
        damage[nodeId] = totalDamage;
    }
}


