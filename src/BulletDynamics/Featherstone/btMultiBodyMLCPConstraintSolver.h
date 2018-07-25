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
#include "LinearMath/btThreads.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"

class btMLCPSolverInterface;
class btMultiBody;

class btMultiBodyMLCPConstraintSolver : public btMultiBodyConstraintSolver
{
protected:
	/// \name MLCP formulation
	/// \{

	/// A matrix in the MLCP formulation
	btMatrixXu m_A;

	/// b vector in the MLCP formulation.
	btVectorXu m_b;

	/// Constraint impulse, which is an output of MLCP solving.
	btVectorXu m_x;

	/// Lower bound of constraint impulse, \c m_x.
	btVectorXu m_lo;

	/// Upper bound of constraint impulse, \c m_x.
	btVectorXu m_hi;

	/// \}

	/// \name Cache Variables for Split Impulse
	/// When using 'split impulse' we solve two separate (M)LCPs
	/// \{

	/// Split impulse Cache vector corresponding to \c m_b.
	btVectorXu m_bSplit;

	/// Split impulse cache vector corresponding to \c m_x.
	btVectorXu m_xSplit;

	/// \}

	/// Indices to normal constraint associated with frictional contact constraint.
	///
	/// This is used by the MLCP solver to update the upper bounds of frictional contact impulse given intermediate normal contact impulse.
	btAlignedObjectArray<int> m_limitDependencies;

	/// Array of all the multibody constraints
	btAlignedObjectArray<btMultiBodySolverConstraint*> m_allConstraintPtrArray;

	/// MLCP solver
	btMLCPSolverInterface* m_solver;

	/// Count of fallbacks of using btSequentialImpulseConstraintSolver, which happens when the MLCP solver fails.
	int m_fallback;

	/// Constructs MLCP terms, which are \c m_A, \c m_b, \c m_lo, and \c m_hi.
	virtual void createMLCPFast(const btContactSolverInfo& infoGlobal);

	/// Solves MLCP and returns the success
	virtual bool solveMLCP(const btContactSolverInfo& infoGlobal);
	// Note: Identical to btMLCPSolver::solveMLCP().

	// Documentation inherited
	btScalar solveGroupCacheFriendlySetup(
		btCollisionObject** bodies,
		int numBodies,
		btPersistentManifold** manifoldPtr,
		int numManifolds,
		btTypedConstraint** constraints,
		int numConstraints,
		const btContactSolverInfo& infoGlobal,
		btIDebugDraw* debugDrawer) BT_OVERRIDE;

	// Documentation inherited
	btScalar solveGroupCacheFriendlyIterations(
		btCollisionObject** bodies ,
		int numBodies,btPersistentManifold** manifoldPtr,
		int numManifolds,btTypedConstraint** constraints,
		int numConstraints,
		const btContactSolverInfo& infoGlobal,
		btIDebugDraw* debugDrawer) BT_OVERRIDE;

public:
	BT_DECLARE_ALIGNED_ALLOCATOR()

	/// Constructor
	///
	/// \param[in] solver MLCP solver. Assumed it's not null.
	explicit btMultiBodyMLCPConstraintSolver(btMLCPSolverInterface* solver);

	/// Destructor
	virtual ~btMultiBodyMLCPConstraintSolver();

	/// Sets MLCP solver
	void setMLCPSolver(btMLCPSolverInterface* solver);

	/// Returns the number of fallbacks of using btSequentialImpulseConstraintSolver, which happens when the MLCP solver fails.
	int getNumFallbacks() const;

	/// Sets the number of fallbacks. This function may be used to reset the number to zero.
	void setNumFallbacks(int num);

	/// Returns the constraint solver type.
	virtual btConstraintSolverType getSolverType() const;
};

#endif  // BT_MULTIBODY_MLCP_CONSTRAINT_SOLVER_H
