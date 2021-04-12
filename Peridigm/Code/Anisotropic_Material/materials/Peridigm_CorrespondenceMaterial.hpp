//! \file Peridigm_CorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_CORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_CORRESPONDENCEMATERIAL_HPP

#include "Peridigm_Material.hpp"
#include "Peridigm_InfluenceFunction.hpp"

namespace PeridigmNS {

  class CorrespondenceMaterial : public Material{
  public:

    //! Constructor.
    CorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~CorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Correspondence Base Class"); }

    //! Returns the density of the material.
    virtual double Density() const { return m_density; }

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const { return m_bulkModulus; }

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const { return m_shearModulus; }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Initialize the material model.
    virtual void initialize(const double dt,
                            const int numOwnedPoints,
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager);

    //! Evaluate the Cauchy stress (pure virtual function, must be implemented by derived correspondence material models).
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     PeridigmNS::DataManager& dataManager) const = 0;

    //! Evaluate the internal force.
    virtual void computeForce(const double dt,
                              const int numOwnedPoints,
                              const int* ownedIDs,
                              const int* neighborhoodList,
                              PeridigmNS::DataManager& dataManager) const;
//////////////////////////////////////////////////////////////////////////////////

//! Evaluate the jacobian.
    virtual void
    computeJacobian(const double dt,
                    const int numOwnedPoints,
                    const int* ownedIDs,
                    const int* neighborhoodList,
                    PeridigmNS::DataManager& dataManager,
                    PeridigmNS::SerialMatrix& jacobian,
                    PeridigmNS::Material::JacobianType jacobianType = PeridigmNS::Material::FULL_MATRIX) const;
    /// \enum JacobianType
    /// \brief Whether to compute the full tangent stiffness matrix or just its block diagonal entries
    ///
    /// The Peridigm Material base class provides a computeJacobian method that all materials inherit
    /// to compute the tangent stiffness matrix. This base class uses a finite difference method (it provides both
    /// forward and centered) to numerically approximate the jacobian. Derived classes may override this method
    /// to compute the jacobian via another approach (for example, automatic differentiation) or may simply
    /// inherit and use the finite difference Jacobian, which will work for all derived mateiral classes.
    ///
    /// The default behavior of this
    /// method is to compute the full tangent stiffness matrix. However, is it occasionally useful to compute
    /// and store only the block diagonal entries -- specifically, the only the entries of the matrix that
    /// describe interactions between different dofs for an individual node. This manifests as a block-diagonal
    /// with each block of dimension 3x3, and one block for each node. This block-diagonal matrix is useful,
    /// for example, as a preconditioner for the full Jacobian matrix. The 3x3 blocks are also the "P" matrices
    /// used by the ComputeStabilityIndex compute class. See that class for more information.
    ///
    /// \note The default behavior is to compute the full tangent stiffness matrix. This enum is useful to only
    /// if you need to efficiently compute only the block diagonal entries of the full tangent stiffness matrix.
    //enum JacobianType { UNDEFINED=0, NONE=1, FULL_MATRIX=2, BLOCK_DIAGONAL=3 };                           
//! Evaluate the jacobian via automatic differentiation.
    virtual void
    computeAutomaticDifferentiationJacobian(const double dt,
                                            const int numOwnedPoints,
                                            const int* ownedIDs,
                                            const int* neighborhoodList,
                                            PeridigmNS::DataManager& dataManager,
                                            PeridigmNS::SerialMatrix& jacobian,
                                            PeridigmNS::Material::JacobianType jacobianType = PeridigmNS::Material::FULL_MATRIX) const;
//! Evaluate the jacobian via finite difference scheme .
    virtual void
    computeJacobianFiniteDifference(const double dt,
                                    const int numOwnedPoints,
                                    const int* ownedIDs,
                                    const int* neighborhoodList,
                                    PeridigmNS::DataManager& dataManager,
                                    PeridigmNS::SerialMatrix& jacobian,
                                    FiniteDifferenceScheme finiteDifferenceScheme,
                                    PeridigmNS::Material::JacobianType jacobianType = PeridigmNS::Material::FULL_MATRIX) const;
  
    
  
    //! Returns a vector of field IDs that need to be synchronized across block boundaries and MPI boundaries after precompute().
  //  virtual std::vector<int> FieldIdsForSynchronizationAfterPrecompute() const {
  //    std::vector<int> fieldIds;
  //    fieldIds.push_back(m_piolaStressTimesInvShapeTensorXId);
  //    fieldIds.push_back(m_piolaStressTimesInvShapeTensorYId);
  //    fieldIds.push_back(m_piolaStressTimesInvShapeTensorZId);
  //    return fieldIds;
  //  }
    //enum FiniteDifferenceScheme { FORWARD_DIFFERENCE=0, CENTRAL_DIFFERENCE=1 };
  protected:

    // material parameters
    double m_bulkModulus;
    double m_shearModulus;
    double m_density;
    bool m_applyAutomaticDifferentiationJacobian;
    double D;
    int    m_stabilizationType;
    double C[6][6];
    bool   m_planeStress, m_planeStrain;
    bool   m_plane = false;
    bool   nonLin = false;
	bool   lin = true;
	bool   avg = false;
    double m_hourglassCoefficient;
    double scal;
    PeridigmNS::InfluenceFunction::functionPointer m_OMEGA;

    // field spec ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_horizonFieldId;
    int m_volumeFieldId;
    int m_modelCoordinatesFieldId;
    int m_coordinatesFieldId;
    int m_velocitiesFieldId;
    int m_hourglassForceDensityFieldId;
    int m_forceDensityFieldId;
    int m_bondDamageFieldId;
    int m_deformationGradientFieldId;
    int m_shapeTensorInverseFieldId;
    int m_cauchyStressFieldId;
    int m_leftStretchTensorFieldId;
    int m_rotationTensorFieldId;
    int m_unrotatedCauchyStressFieldId;
    int m_piolaStressTimesInvShapeTensorXId;
    int m_piolaStressTimesInvShapeTensorYId;
    int m_piolaStressTimesInvShapeTensorZId;
    int m_unrotatedRateOfDeformationFieldId;
    int m_detachedNodesFieldId;
    int m_partialStressFieldId;
    int m_hourglassStiffId;
    int m_netdamageFieldId;
    int m_type;
    int m_modelAnglesId;
    int m_incremental;
    int m_deformationGradientNonIncFieldId;
    bool m_tension;
  };
}

#endif // PERIDIGM_CORRESPONDENCEMATERIAL_HPP
