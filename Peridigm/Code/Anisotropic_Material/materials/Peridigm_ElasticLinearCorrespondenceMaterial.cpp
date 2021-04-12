/*! \file Peridigm_ElasticLinearCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticLinearCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticLinearCorrespondenceMaterial::ElasticLinearCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_cauchyStressFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_modelAnglesId(-1)
{
  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
  m_density = params.get<double>("Density");
  if (params.isParameter("Plane Strain")){
    m_planeStrain = params.get<bool>("Plane Strain");
    }
  if (params.isParameter("Plane Stress")){
    m_planeStress = params.get<bool>("Plane Stress");
    }
  if (m_planeStrain==true)m_type=1;
  if (m_planeStress==true)m_type=2;
  m_incremental = false;
    if (params.isParameter("Incremental")){
        m_incremental = params.get<bool>("Incremental");
    }
     
 
  double C11=0.0, C44=0.0, C55=0.0, C66=0.0, C12=0.0, C13=0.0, C14=0.0, C15=0.0, C16=0.0, C22=0.0, C33=0.0;
  double  C23=0.0, C24=0.0, C25=0.0, C26=0.0, C34=0.0, C35=0.0, C36=0.0, C45=0.0, C46=0.0, C56=0.0;
  
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
       //std::cout<<"??"<<std::endl;
   }
  }
  else{
       m_bulkModulus = calculateBulkModulus(params);
       m_shearModulus = calculateShearModulus(params);
  
       C11 = 2*m_shearModulus + m_bulkModulus;
       C44 = m_shearModulus;
       C55 = m_shearModulus;
       C66 = m_shearModulus;
       C12 = C11 - 2*C55;
       C13 = C12;
       C14 = 0.0;
       C15 = 0.0;
       C16 = 0.0;
       C22 = 2*m_shearModulus + m_bulkModulus;
       C33 = 2*m_shearModulus + m_bulkModulus;
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
  if (m_plane==false){
   C[0][0] = C11;C[0][1] = C12;C[0][2]= C13; C[0][3] = C14; C[0][4] = C15; C[0][5]= C16;
   C[1][0] = C12;C[1][1] = C22;C[1][2]= C23; C[1][3] = C24; C[1][4] = C25; C[1][5]= C26;
   C[2][0] = C13;C[2][1] = C23;C[2][2]= C33; C[2][3] = C34; C[2][4] = C35; C[2][5]= C36;
   C[3][0] = C14;C[3][1] = C24;C[3][2]= C34; C[3][3] = C44; C[3][4] = C45; C[3][5]= C46;
   C[4][0] = C15;C[4][1] = C25;C[4][2]= C35; C[4][3] = C45; C[4][4] = C55; C[4][5]= C56;
   C[5][0] = C16;C[5][1] = C26;C[5][2]= C36; C[5][3] = C46; C[5][4] = C56; C[5][5]= C66;
  }
  // tbd in future
  if (m_planeStress==true){
  //    //only transversal isotropic in the moment
   C[0][0] = C11-C13*C13/C22;C[0][1] = C12-C13*C23/C22;C[0][2] = 0.0; C[0][3] = 0.0; C[0][4] = 0.0; C[0][5] = 0.0;
   C[1][0] = C12-C13*C23/C22;C[1][1] = C22-C13*C23/C22;C[1][2] = 0.0; C[1][3] = 0.0; C[1][4] = 0.0; C[1][5] = 0.0;
   C[2][0] = 0.0;            C[2][1] = 0.0;            C[2][2] = 0.0; C[2][3] = 0.0; C[2][4] = 0.0; C[2][5] = 0.0;
   C[3][0] = 0.0;            C[3][1] = 0.0;            C[3][2] = 0.0; C[3][3] = 0.0; C[3][4] = 0.0; C[3][5] = 0.0;
   C[4][0] = 0.0;            C[4][1] = 0.0;            C[4][2] = 0.0; C[4][3] = 0.0; C[4][4] = 0.0; C[4][5] = 0.0;
   C[5][0] = 0.0;            C[5][1] = 0.0;            C[5][2] = 0.0; C[5][3] = 0.0; C[5][4] = 0.0; C[5][5] = C66;
   
  
  }
  // not correct for plane stress!!!
  if (m_plane!=0){
   C[0][0] = C11;C[0][1] = C12;C[0][2] = 0.0;C[0][3] = 0.0;C[0][4] = 0.0;C[0][5] = C16;
   C[1][0] = C12;C[1][1] = C22;C[1][2] = 0.0;C[1][3] = 0.0;C[1][4] = 0.0;C[1][5] = C26;
   C[2][0] = 0.0;C[2][1] = 0.0;C[2][2] = 0.0;C[2][3] = 0.0;C[2][4] = 0.0;C[2][5] = 0.0;
   C[3][0] = 0.0;C[3][1] = 0.0;C[3][2] = 0.0;C[3][3] = 0.0;C[3][4] = 0.0;C[3][5] = 0.0;
   C[4][0] = 0.0;C[4][1] = 0.0;C[4][2] = 0.0;C[4][3] = 0.0;C[4][4] = 0.0;C[4][5] = 0.0;
   C[5][0] = C16;C[5][1] = C26;C[5][2] = 0.0;C[5][3] = 0.0;C[5][4] = 0.0;C[5][5] = C66;
  }
  
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();

  m_cauchyStressFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_deformationGradientFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
  m_modelAnglesId                       = fieldManager.getFieldId(PeridigmField::NODE   , PeridigmField::VECTOR, PeridigmField::CONSTANT     , "Local_Angles");
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
 


}

PeridigmNS::ElasticLinearCorrespondenceMaterial::~ElasticLinearCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticLinearCorrespondenceMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager)
{
      PeridigmNS::CorrespondenceMaterial::initialize(dt,
                                                      numOwnedPoints,
                                                      ownedIDs,
                                                      neighborhoodList,
                                                      dataManager);

                             
}

void
PeridigmNS::ElasticLinearCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               PeridigmNS::DataManager& dataManager) const
{

  

  double *CauchyStress, *CauchyStressNP1, *defGrad, *angles;
  // have to be checked if the additional effort is useful or not
  // deactivated for tests with implicit solver
  bool incremental = false;
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&CauchyStress);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);

  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&defGrad);
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);

  CORRESPONDENCE::updateElasticCauchyStressSmallDef(defGrad, 
                                            CauchyStress,
                                            CauchyStressNP1,
                                            numOwnedPoints,
                                            C,
                                            angles,
                                            m_type,
                                            dt,
                                            incremental);
                                            
    if (m_incremental == true){
        std::cout<<"please choose elastic correspondence. incremental not implemented yet for anisotropic material"<<std::endl;
    }                                        
                                            
                                           
}
