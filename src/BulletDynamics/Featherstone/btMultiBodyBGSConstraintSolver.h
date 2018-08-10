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

#ifndef BT_MULTIBODY_BGS_CONSTRAINT_SOLVER_H
#define BT_MULTIBODY_BGS_CONSTRAINT_SOLVER_H

#include "LinearMath/btMatrixX.h"
#include "LinearMath/btThreads.h"
#include "BulletDynamics/MLCPSolvers/btMLCP.h"
#include "BulletDynamics/MLCPSolvers/btSolveProjectedGaussSeidel.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyMLCP.h"

class btMLCPSolverInterface;
class btMultiBody;

class btMultiBodyBGSConstraintSolver : public btMultiBodyConstraintSolver
{
protected:
	/// Array of MLCP blocks for rigid bodies.
	btAlignedObjectArray<btMLCP> m_mlcpArray;
	// Note: If we know the size of MLCP in compile time, we could use a fixed size MLCP struct, which may don't
	// require memory allocations in the simulation loops.

	/// Array of MLCP blocks for multibodies.
	btAlignedObjectArray<btMultiBodyMLCP> m_multiBodyMlcpArray;
	// Note: If we know the size of MLCP in compile time, we could use a fixed size MLCP struct, which may don't
	// require memory allocations in the simulation loops.

	/// Default MLCP solver
	btMLCPSolverInterface* m_defaultSolver;

	/// PGS MLCP solver to be used when \c m_defaultSolver is failed.
	btSolveProjectedGaussSeidel m_pgsSolver;

	/// Count of fallbacks of using PGS LCP solver (\c m_pgsSolver), which happens when the default MLCP solver (\c m_defaultSolver) fails.
	int m_fallback;

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

	/// Constructs a MLCP for the constraints associated with a contact point such as normal contact constraints,
	/// frictional constraints, and torsional constraints. The constructed MLCP block is stored in \c m_mlcpArray.
	///
	/// \param[in] index Index to the MLCP block.
	/// \param[in] numFrictionPerContact Number of friction constraints per contact.
	/// \param[in] infoGlobal Global configurations for contact solver.
	void setupContactConstraintMLCPBlockRigidBody(
		int index,
		int numFrictionPerContact,
		const btContactSolverInfo& infoGlobal);

	void setupContactConstraintMLCPBlockMultiBody(
		int index,
		int numFrictionPerContact,
		const btContactSolverInfo& infoGlobal);

	// Documentation inherited.
	btScalar solveSingleIteration(
		int iteration,
		btCollisionObject** bodies,
		int numBodies,
		btPersistentManifold** manifoldPtr,
		int numManifolds,
		btTypedConstraint** constraints,
		int numConstraints,
		const btContactSolverInfo& infoGlobal,
		btIDebugDraw* debugDrawer) BT_OVERRIDE;

	/// Solves a MLCP block, which are stored in \c m_mlcpArray.
	///
	/// \param[in] index Index to the MLCP block.
	/// \param[in] infoGlobal Global configurations for contact solver.
	btScalar solveMLCPBlockRigidBody(int index, const btContactSolverInfo& infoGlobal);

	btScalar solveMLCPBlockMultiBody(int index, const btContactSolverInfo& infoGlobal);

public:
	BT_DECLARE_ALIGNED_ALLOCATOR()

	/// Constructor
	///
	/// \param[in] solver MLCP solver. Assumed it's not null.
	explicit btMultiBodyBGSConstraintSolver(btMLCPSolverInterface* solver);

	/// Destructor
	virtual ~btMultiBodyBGSConstraintSolver();

	/// Sets MLCP solver. Assumed it's not null.
	void setMLCPSolver(btMLCPSolverInterface* solver);

	/// Returns the number of fallbacks of using btSequentialImpulseConstraintSolver, which happens when the MLCP
	/// solver fails.
	int getNumFallbacks() const;

	/// Sets the number of fallbacks. This function may be used to reset the number to zero.
	void setNumFallbacks(int num);

	/// Returns the constraint solver type.
	virtual btConstraintSolverType getSolverType() const;
};

#endif  // BT_MULTIBODY_BGS_CONSTRAINT_SOLVER_H
