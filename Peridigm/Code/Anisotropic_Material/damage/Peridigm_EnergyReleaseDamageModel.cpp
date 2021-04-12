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
// https://github.com/PeriDoX/
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
m_planeStrain(false),
m_planeStress(false),
m_onlyTension(false),
m_rot(false),
m_Thickness(-1),
m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()) {
	m_rot = false;
    if (params.isParameter("Rot Sym")) {
		m_rot = true;	
		maxRad = params.get<double>("Rot Sym");
	}
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
        if (params.isParameter("Energy Criterion")) {
            if (params.get<bool>("Energy Criterion") == true)
                m_type = 1;
            }
            
        if (params.isParameter("Power Law")) {
            if (params.get<bool>("Power Law") == true)
                m_type = 2;
            }
            
        if (params.isParameter("Separated")) {
            if (params.get<bool>("Separated") == true)
                m_type = 3;
            }

    }
    if(params.isParameter("Plane Stress")){
        m_planeStress = params.get<bool>("Plane Stress");
        m_Thickness = params.get<double>("Thickness");
    }
    if(params.isParameter("Plane Strain")){
        m_planeStrain = params.get<bool>("Plane Strain");
        m_Thickness = params.get<double>("Thickness");
    }
    m_pi = 3.14159;
    m_onlyTension = false;
    if(params.isParameter("Only Tension")){
        m_onlyTension = params.get<bool>("Only Tension");
       
    }

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
    if(m_applyThermalStrains)
        m_deltaTemperatureFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Temperature_Change");
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
	dataManager.getData(m_damageModelFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
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
    if(m_applyThermalStrains)
        dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
    double m_alpha = 0;

    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
    //////////////////////////////////////////////////////////////
    // transfer of data is done in ComputeDilation --> PeridigmElastic.cpp
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
    double totalDamage, totalVol;
    double alphaP1 = 0.0, alphaP2 = 0.0, gammaP1 = 0.0, gammaP2 = 0.0, kappaP1 = 0.0, kappaP2 = 0.0;
    double nodeInitialX[3], nodeCurrentX[3], relativeExtension(0.0);
    double bondEnergyIsotropic(0.0), bondEnergyDeviatoric(0.0);
    double omegaP1, omegaP2;
    double critDev, critIso, critComp;
    double degradationFactor = 1.0; // Optional parameter if bond should be degradated and not fully destroyed instantaneously 
    double avgHorizon, quadhorizon;
    double eiP1, eiP2, tiP1, tiP2;
    double edP1, edP2, tdP1, tdP2;
    double dilatationP1, dilatationP2;
    double BulkModP1, BulkModP2;
    double ShearModP1, ShearModP2;
    double bondEnergy;
    
    //double thickness;
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
    //temperatur ist in peridigm_ElasticMaterial gegeben
    for (iID = 0; iID < numOwnedPoints; ++iID) {
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        
        nodeId = ownedIDs[iID];
        nodeInitialX[0] = x[nodeId*3];
        nodeInitialX[1] = x[nodeId*3+1];
        nodeInitialX[2] = x[nodeId*3+2];
        nodeCurrentX[0] = y[nodeId*3];
        nodeCurrentX[1] = y[nodeId*3+1];
        nodeCurrentX[2] = y[nodeId*3+2];
 
        dilatationP1 = damageModel[3*nodeId];
        BulkModP1    = damageModel[3*nodeId+1];     // weightedVolume is already included in the variable --> Perdigm_Elastic.cpp;
        ShearModP1   = damageModel[3*nodeId+2]; // weightedVolume is already included in --> Perdigm_Elastic.cpp;
        
        if (m_planeStress==false and m_planeStrain==false){
          //  c = 3 * K * (*theta) / weightedVol;
            // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
            kappaP1 = 9.0 * BulkModP1;
            gammaP1 = 1.0;
            alphaP1 = 15.0 * ShearModP1;
        }
        if (m_planeStress==true and m_planeStrain==false){
            alphaP1 = 8.0 * ShearModP1;
            kappaP1 = 3.0 * BulkModP1;
            gammaP1 = 4.0 * ShearModP1 / (3.0 * BulkModP1 + 4.0 * ShearModP1);
            alphaP1 = 8.0 * ShearModP1;
        }
        if (m_planeStrain==true and m_planeStress==false){
           // c = (12.0*K-4.0*MU) / 9.0 * (*theta) / weightedVol;
            // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
            //kappa = 3*(12.0*K-4.0*MU) / 18.0 * (*theta) / weightedVol;
            
            gammaP1 = 2.0/3.0;
            kappaP1 = 6.0*BulkModP1 - 2.0*ShearModP1;
            alphaP1 = 8.0*ShearModP1;
        }
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
            if(deltaTemperature)
              eP -= 0.5 * m_alpha*(deltaTemperature[nodeId] + deltaTemperature[neighborID])*zeta;
            
            
            // the average horizon is taken, if multiple horizons are used
            dilatationP2 = damageModel[3*neighborID];
            BulkModP2    = damageModel[3*neighborID+1]; // weightedVolume is already included in the variable --> Perdigm_Elastic.cpp;
            ShearModP2   = damageModel[3*neighborID+2]; // weightedVolume is already included in --> Perdigm_Elastic.cpp;
			
            
            
            
            if (m_planeStress){
              //  c = 4.0*K*MU/(3.0*K+4.0*MU) * (*theta) / weightedVol;
                // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
                kappaP2 = 3.0 * BulkModP2;
                gammaP2 = 4.0 * ShearModP2 / (3.0 * BulkModP2 + 4.0 * ShearModP2);
                alphaP2 = 8.0 * ShearModP2;
            }
            if (m_planeStrain){
               // c = (12.0*K-4.0*MU) / 9.0 * (*theta) / weightedVol;
                // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
                //kappa = 3*(12.0*K-4.0*MU) / 18.0 * (*theta) / weightedVol;
                
                gammaP2 = 2.0/3.0;
                kappaP2 =(12.0*BulkModP2 - 4.0*ShearModP2) / 3.0 / gammaP2;
                alphaP2 = 8.0 * ShearModP2;
                
            }
            if (m_planeStress==false and m_planeStrain==false){
              //  c = 3 * K * (*theta) / weightedVol;
                // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
                kappaP2 = 9.0 * BulkModP2;
                gammaP2 = 1.0;
                alphaP2 = 15.0 * ShearModP2 ;
            }

            omegaP1 = MATERIAL_EVALUATION::scalarInfluenceFunction(zeta, horizon[nodeId]); 
            omegaP2 = MATERIAL_EVALUATION::scalarInfluenceFunction(zeta, horizon[neighborID]); 
            avgHorizon = 0.5*(horizon[nodeId]+horizon[neighborID]);
            
            if(m_planeStrain or m_planeStress){ // m_Thickness is in the volume included (tbd to guarantee this)
            //factor = 16.0 / 5.0 / weightedVol * sqrt(M_PI / 3.0 * pow(horizon,5)); // Simplified with Mathematica
                quadhorizon =  3.0 /(m_Thickness * pow(horizon[neighborID],3) );
				if (m_rot){
					double rad = nodeInitialX[0]*nodeInitialX[0]/maxRad/maxRad;
					//if (rad == 0) rad = 1e-10;
					quadhorizon = quadhorizon*(1+rad);
				}
            }
            else {
                //yieldValue = 25.0 * yieldStress * yieldStress / 8 / M_PI / pow(horizon,5);
                //factor = sqrt(75.0 / 4 / M_PI / pow(horizon,5));        // equation (51) MitchelJA_2011; 1/2 is not needed 
                //factor = 15.0 / weightedVol * sqrt(4.0 / 75.0 * M_PI * pow(horizon,5));        // equation (51) MitchelJA_2011; 1/2 is not needed 
                
                quadhorizon =  8.0 /( m_pi * pow(avgHorizon, 4) );
            }
            


            eiP1 = gammaP1 * dilatationP1 * zeta / 3.0;
            eiP2 = gammaP2 * dilatationP2 * zeta / 3.0;
            tiP1 = kappaP1*eiP1;
            tiP2 = kappaP2*eiP2;
            
			
            edP1 = eP - eiP1;
            edP2 = eP - eiP2;
            tdP1 = omegaP1*alphaP1*edP1;
            tdP2 = omegaP2*alphaP2*edP2;
 
            bondEnergyIsotropic  = (1.0 - bondDamageNP1[bondIndex])*(sqrt(tiP1*eiP1*tiP1*eiP1) + sqrt(tiP2*eiP2*tiP2*eiP2)); 
            bondEnergyDeviatoric = (1.0 - bondDamageNP1[bondIndex])*(sqrt(tdP1*edP1*tdP1*edP1) + sqrt(tdP2*edP2*tdP2*edP2));

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
                bondEnergy = abs(abs(tdP1+tiP1) + abs(tdP2 + tiP2))*eP;
                
                //bondEnergy = abs(damageModel[3*nodeId]-damageModel[3*neighborID])*zeta / *m;
                
                critIso = bondEnergy/(criticalEnergyTension*quadhorizon);
                
                if (m_criticalEnergyTension > 0.0 && critIso > 1.0 ) {
					//std::cout<<BulkModP1<< " "<< weightedVolume[nodeId]<<" "<< ShearModP1<< " "<< ShearModP2<<" "<< BulkModP2<<std::endl;
                    trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                } 
            }
            
            if (m_type == 2){// Power Law
    
                if (m_criticalEnergyTension > 0.0 &&critIso*critIso + critDev*critDev + critComp*critComp > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
            }
            //std::cout<< m_type<< " "<< std::endl;
            if (m_type == 3){ // Separated
                if (m_criticalEnergyTension > 0.0  &&critIso > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
                if (criticalEnergyShear > 0.0 && critDev > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                   
                } 
                if (m_criticalEnergyCompression > 0.0  && critComp > 1.0) {
                   trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
            }
            if (m_onlyTension == true){
                if (relativeExtension<0)trialDamage = 0;
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
        totalVol = 0.0;
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            
            neighborID = neighborhoodList[neighborhoodListIndex++];
            // must be zero to avoid synchronization errors
            damageModel[3*neighborID] = 0.0;
            damageModel[3*neighborID+1] = 0.0;
            damageModel[3*neighborID+2] = 0.0;
            
            totalDamage += bondDamageNP1[bondIndex]*weightedVolume[nodeId];
            totalVol += weightedVolume[nodeId];
            bondIndex += 1;
        }
        if (numNeighbors > 0)
            totalDamage /= (totalVol);
        else
            totalDamage = 0.0;
        damage[nodeId] = totalDamage;
    }
}


