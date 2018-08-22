/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2018 Google Inc. http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use
of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim
that you wrote the original software. If you use this software in a product, an
acknowledgment in the product documentation would be appreciated but is not
required.
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef BT_BLOCK_GAUSS_SEIDEL_CONSTRAINT_SOLVER_H
#define BT_BLOCK_GAUSS_SEIDEL_CONSTRAINT_SOLVER_H

#include "BulletDynamics/ConstraintSolver/btConstraintSolver.h"
#include "BulletDynamics/ConstraintSolver/btContactSolverInfo.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btThreads.h"

class btConstraintBlockGenerator;

struct btConstraintBlock
{
	btSolverConstraintInput m_input;
	btConstraintSolver* m_solver;

	/// Constructor
	btConstraintBlock(const btSolverConstraintInput& input, btConstraintSolver* solver);

	btScalar solve(const btContactSolverInfo& infoGlobal, class btIDebugDraw* debugDrawer, btDispatcher* dispatcher);
};

/// Implementation of the blocked Gauss-Seidel (BGS) constraint solver.
///
/// Loosely speaking this class is a mix of btSequentialImpulseConstraintSolver
/// and btMLCPSolver. Specifically, BSG
/// forms many small MLCPs of subsets of the constraints and solves each MLCP
/// using solvers suitable for small-sized
/// problems (e.g., Dantzig), and then applies the Gauss-Seidel splitting to the
/// blocked MLCPs unlike the regular
/// Gauss-Seidel applies it to every single constraint. A block is usually
/// defined as all constraints associated with a
/// single contact point or any set of constraints better to be solved together
/// at once.
///
/// The expected performance is somewhere between
/// btSequentialImpulseConstraintSolver and btMLCPSolver in terms of
/// speed and accuracy.
class btBlockGaussSeidelConstraintSolver : public btConstraintSolver
{
protected:
	btConstraintBlockGenerator* m_splittingPolicy;

	btAlignedObjectArray<btConstraintBlock> m_blocks;

	int m_maxOverrideNumSolverIterations;

	btScalar m_squaredResidual;

	virtual btScalar solveGroupCacheFriendlySetup(
		btCollisionObject** bodies, int numBodies,
		btPersistentManifold** manifoldPtr, int numManifolds,
		btTypedConstraint** constraints, int numConstraints,
		const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer,
		btDispatcher* dispatcher);

	virtual btScalar
	solveGroupCacheFriendlyIterations(const btContactSolverInfo& infoGlobal,
									  btIDebugDraw* debugDrawer,
									  btDispatcher* dispatcher);

	virtual btScalar solveSingleIteration(int iteration,
										  const btContactSolverInfo& infoGlobal,
										  btIDebugDraw* debugDrawer,
										  btDispatcher* dispatcher);

	virtual btScalar
	solveGroupCacheFriendlyFinish(btCollisionObject** bodies, int numBodies,
								  const btContactSolverInfo& infoGlobal);

public:
	btBlockGaussSeidelConstraintSolver(btConstraintBlockGenerator* splittingPolicy = 0);

	~btBlockGaussSeidelConstraintSolver() BT_OVERRIDE;

	void setSplittingPolicy(btConstraintBlockGenerator* splittingPolity);

	btConstraintBlockGenerator* getSplittingPolicy();

	const btConstraintBlockGenerator* getSplittingPolicy() const;

	// Documentation inherited
	btScalar solveGroup(btCollisionObject** bodies, int numBodies,
						btPersistentManifold** manifoldPtr, int numManifolds,
						btTypedConstraint** constraints, int numConstraints,
						const btContactSolverInfo& infoGlobal,
						class btIDebugDraw* debugDrawer,
						btDispatcher* dispatcher) BT_OVERRIDE;

	// Documentation inherited
	void reset() BT_OVERRIDE;

	btConstraintSolverType getSolverType() const BT_OVERRIDE;
};

#endif  // BT_BLOCK_GAUSS_SEIDEL_CONSTRAINT_SOLVER_H
