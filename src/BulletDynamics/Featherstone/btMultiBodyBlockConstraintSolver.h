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

#ifndef BT_MULTIBODY_BLOCK_CONSTRAINT_SOLVER_H
#define BT_MULTIBODY_BLOCK_CONSTRAINT_SOLVER_H

#include "btMultiBodyConstraintSolver.h"

#include "LinearMath/btThreads.h"

struct btBlockConstraintSolverConfig
{
	int m_solverType;  //SI or MLCP Dantzig
	//to be decided: full or subset of

	btContactSolverInfo m_info;
};

struct btConstraintBlock
{
	/// \{ \name Rigid Body Data

	btCollisionObject** m_bodies;
	int m_numBodies;
	btPersistentManifold** m_manifold;
	int m_numManifolds;
	btTypedConstraint** m_constraints;
	int m_numConstraints;
	btMultiBodyConstraint** m_multiBodyConstraints;
	int m_numMultiBodyConstraints;

	/// Pointer to the block constraint solver's body pool, which will be shared
	/// by all the constraint blocks.
	btAlignedObjectArray<btSolverBody>* m_tmpSolverBodyPool;

	btConstraintArray m_tmpSolverContactConstraintPool;
	btConstraintArray m_tmpSolverNonContactConstraintPool;
	btConstraintArray m_tmpSolverContactFrictionConstraintPool;
	btConstraintArray m_tmpSolverContactRollingFrictionConstraintPool;

	/// \}

	/// \{ \name Multibody Data

	btMultiBodyConstraintArray m_multiBodyNonContactConstraints;

	btMultiBodyConstraintArray m_multiBodyNormalContactConstraints;
	btMultiBodyConstraintArray m_multiBodyFrictionContactConstraints;
	btMultiBodyConstraintArray m_multiBodyTorsionalFrictionContactConstraints;
	// TODO(JS): Change the array names to be more consistent for both of rigid body and multibody

	/// Pointer to the block constraint solver's multi body Jacobian data, which
	/// will be shared by all the constraint blocks.
	btMultiBodyJacobianData* m_data;

	btMultiBodyConstraint** m_tmpMultiBodyConstraints;
	int m_tmpNumMultiBodyConstraints;

	/// \}

	btMultiBodyConstraintSolver* m_solver;

	int m_constraintConfigId;

	btConstraintBlock();
};

class btBlockSplittingPolicy
{
protected:
	virtual ~btBlockSplittingPolicy();

public:
	/// Splits a set of constraints into multiple subsets.
	///
	/// \param[in] blockInput
	/// \param[in] availableConfigs
	/// \param[in,out] blocksOutput
	virtual void split(
		const btConstraintBlock& blockInput,
		const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs,
		btAlignedObjectArray<btConstraintBlock>& blocksOutput) = 0;
};

class btMultiBodyBlockConstraintSolver : public btMultiBodyConstraintSolver
{
protected:
	/// Splitting policy. Assumed not a null.
	btBlockSplittingPolicy* m_splittingPolicy;

	/// Array of constraint configurations for constraint blocks.
	btAlignedObjectArray<btBlockConstraintSolverConfig> m_configs;

	/// Array of constraint blocks.
	btAlignedObjectArray<btConstraintBlock> m_blocks;

public:
	/// Constructor
	btMultiBodyBlockConstraintSolver();

	/// Destructor
	virtual ~btMultiBodyBlockConstraintSolver();

protected:
	// Documentation inherited.
	void solveMultiBodyGroup(btCollisionObject** bodies,
							 int numBodies,
							 btPersistentManifold** manifold,
							 int numManifolds,
							 btTypedConstraint** constraints,
							 int numConstraints,
							 btMultiBodyConstraint** multiBodyConstraints,
							 int numMultiBodyConstraints,
							 const btContactSolverInfo& info,
							 btIDebugDraw* debugDrawer,
							 btDispatcher* dispatcher) BT_OVERRIDE;

	/// Sets the splitting policy.
	virtual void setSplittingPolicy(btBlockSplittingPolicy* policy);

	/// Adds a constraint block configuration and returns the total number of configurations added to this solver.
	virtual int addConfig(btBlockConstraintSolverConfig& config);

	/// Returns the number of configurations added to this solver.
	virtual int getNumConfigs() const;

	/// Removes an configuration at \c configIndex
	///
	/// \param[in] configIndex The configuration indext in the range of [0, numConfigs). Passing out of the range is an
	/// undefined behavior.
	virtual void removeConfig(int configIndex);  //in range 0..numConfigs
};

#endif  //BT_MULTIBODY_BLOCK_CONSTRAINT_SOLVER_H
