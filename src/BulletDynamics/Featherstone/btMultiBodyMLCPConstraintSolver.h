/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2018 Google Inc. http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef BT_MULTIBODY_MLCP_CONSTRAINT_SOLVER_H
#define BT_MULTIBODY_MLCP_CONSTRAINT_SOLVER_H

#include "LinearMath/btMatrixX.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "BulletDynamics/MLCPSolvers/btMLCPSolverInterface.h"

#define DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS

class btMultiBody;

// TODO(JS): btRigidBody vs btRigidBody isn't handled for now
class btMultiBodyMLCPConstraintSolver : public btMultiBodyConstraintSolver
{
protected:
	/// \name Boxed LCP formulation
	/// \{

	// TODO(JS): Consider using Eigen, which is suggested by Erwin
	btMatrixXu m_A;
	btVectorXu m_b;
	btVectorXu m_x;
	btVectorXu m_lo;
	btVectorXu m_hi;

	/// \}

	/// \name Cache Variables for Split Impulse
	/// When using 'split impulse' we solve two separate (M)LCPs
	/// \{

	btVectorXu m_bSplit;
	btVectorXu m_xSplit;
	btVectorXu m_bSplit1;
	btVectorXu m_xSplit2;

	/// \}

	btAlignedObjectArray<int> m_limitDependencies;
	btAlignedObjectArray<btSolverConstraint*> m_allConstraintPtrArray;
	btMLCPSolverInterface* m_solver;
	int m_fallback;

	/// \name Scratch Variables
	/// The following scratch variables are not stateful -- contents are cleared prior to each use.
	/// They are only cached here to avoid extra memory allocations and deallocations and to ensure
	/// that multiple instances of the solver can be run in parallel.
	/// \{

	btMatrixXu m_scratchJ3;
	btMatrixXu m_scratchJInvM3;
	btAlignedObjectArray<int> m_scratchOfs;
	btMatrixXu m_scratchMInv;
	btMatrixXu m_scratchJ;
	btMatrixXu m_scratchJTranspose;
	btMatrixXu m_scratchTmp;

	/// \}

	/// Creates Mixed LCP
	virtual void createMLCPFast(const btContactSolverInfo& infoGlobal);

	//return true is it solves the problem successfully
	virtual bool solveMLCP(const btContactSolverInfo& infoGlobal);

	// Documentation inherited
	virtual btScalar solveGroupCacheFriendlySetup(
		btCollisionObject** bodies,
		int numBodies,
		btPersistentManifold** manifoldPtr,
		int numManifolds,
		btTypedConstraint** constraints,
		int numConstraints,
		const btContactSolverInfo& infoGlobal,
		btIDebugDraw* debugDrawer);

	// Documentation inherited
	virtual btScalar solveGroupCacheFriendlyIterations(
		btCollisionObject** bodies ,
		int numBodies,btPersistentManifold** manifoldPtr,
		int numManifolds,btTypedConstraint** constraints,
		int numConstraints,
		const btContactSolverInfo& infoGlobal,
		btIDebugDraw* debugDrawer);

	// Documentation inherited
	void convertContacts(btPersistentManifold** manifoldPtr, int numManifolds,
						 const btContactSolverInfo& infoGlobal) override;

	void convertMultiBodyContact(btPersistentManifold* manifold, const btContactSolverInfo& infoGlobal);

	// Setup m_data here
	void setupMultiBodyContactConstraint(btMultiBodySolverConstraint& solverConstraint, const btVector3& contactNormal,
										 btManifoldPoint& cp, const btContactSolverInfo& infoGlobal,
										 btScalar& relaxation, bool isFriction, btScalar desiredVelocity = 0,
										 btScalar cfmSlip = 0);
	// TODO(JS): Remove jacDiagABInv

	void applyDeltaVee(btScalar* deltaV, btScalar impulse, int velocityIndex, int ndof);

public:
	BT_DECLARE_ALIGNED_ALLOCATOR()

	btMultiBodyMLCPConstraintSolver(btMLCPSolverInterface* solver);

	virtual ~btMultiBodyMLCPConstraintSolver();

	void setMLCPSolver(btMLCPSolverInterface* solver);

	int getNumFallbacks() const;
	void setNumFallbacks(int num);

	virtual btConstraintSolverType	getSolverType() const;
};

#endif  // BT_MULTIBODY_MLCP_CONSTRAINT_SOLVER_H
