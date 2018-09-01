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

#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"

struct btBlockConstraintSolverConfig
{
	int m_solverType;  //SI or MLCP Dantzig
	//to be decided: full or subset of

	btContactSolverInfo m_info;
};

struct btConstraintBlock
{
	// TODO(JS): These are not used actually.
	btCollisionObject** m_bodies;
	int m_numBodies;
	btPersistentManifold** m_manifold;
	int m_numManifolds;

	/// \{ \name Rigid Body Data

	/// Rigid body (joint) constraints. This is shared by all the blocks.
//	btTypedConstraint** m_constraints;

	/// Number of rigid body (joint) constraints. This is shared by all the
	/// blocks.
//	int m_numConstraints;

	/// Pointer to the block constraint solver's body pool, which is shared by
	/// all the constraint blocks.
//	btAlignedObjectArray<btSolverBody>* m_solverBodyPool;

	/// Array of non-contact constraints
//	btConstraintArray m_nonContactConstraints;

	/// Array of normal contact constraints
//	btConstraintArray m_normalContactConstraints;

	/// Array of friction contact constraints
//	btConstraintArray m_frictionContactConstraints;

	/// Array of rolling friction contact constraints
//	btConstraintArray m_rollingFrictionContactConstraints;

	/// \}

	/// \{ \name Multibody Data

	/// Multibody (joint) constraints. This is shared by all the blocks.
//	btMultiBodyConstraint** m_multiBodyConstraints;

	/// Number of multibody (joint) constraints. This is shared by all the
	/// blocks.
//	int m_numMultiBodyConstraints;

	/// Array of multibody non-contact constraints
//	btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyNonContactConstraints;

	/// Array of multibody normal contact constraints
//	btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyNormalContactConstraints;

	/// Array of multibody friction contact constraints
//	btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyFrictionContactConstraints;

	/// Array of multibody rolling friction contact constraints
//	btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyTorsionalFrictionContactConstraints;

	/// \}

	btMultiBodyConstraintSolver::btMultiBodyConstraints m_multiBodyConstraintSet;

	/// Constraint solver
	btMultiBodyConstraintSolver* m_solver;

	/// Index to constraint solver configuration
	int m_constraintConfigId;

	/// Default constructor
	btConstraintBlock();

	/// Constructor
	btConstraintBlock(
		btTypedConstraint** m_constraints,
		int m_numConstraints,
		btAlignedObjectArray<btSolverBody>* m_solverBodyPool,
		btConstraintArray& m_nonContactConstraints,
		btConstraintArray& m_normalContactConstraints,
		btConstraintArray& m_frictionContactConstraints,
		btConstraintArray& m_rollingFrictionContactConstraints,
		btMultiBodyConstraint** m_multiBodyConstraints,
		int m_numMultiBodyConstraints,
		btAlignedObjectArray<btMultiBodySolverConstraint>& m_multiBodyNonContactConstraints,
		btAlignedObjectArray<btMultiBodySolverConstraint>& m_multiBodyNormalContactConstraints,
		btAlignedObjectArray<btMultiBodySolverConstraint>& m_multiBodyFrictionContactConstraints,
		btAlignedObjectArray<btMultiBodySolverConstraint>& m_multiBodyTorsionalFrictionContactConstraints,
		btMultiBodyJacobianData* m_data);
};

class btBlockSplittingPolicy
{
protected:
	void copyContactConstraint(const btConstraintBlock& src, int index, btConstraintBlock& dest);

public:
	/// Destructor
	virtual ~btBlockSplittingPolicy();

	/// Splits a set of constraints into multiple subsets.
	///
	/// \param[in] blockInput
	/// \param[in] availableConfigs
	/// \param[in,out] blocksOutput The splitted blocks. This function adds blocks without clearning the array
	/// beforehand. Clearning the array is the caller's responsibility.
	virtual void split(const btConstraintBlock& blockInput, const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs, btAlignedObjectArray<btConstraintBlock>& blocksOutput) = 0;
};

class btSingleBlockSplittingPolicy : public btBlockSplittingPolicy
{
protected:
	btMultiBodyConstraintSolver* m_solver;

public:
	/// Constructor
	btSingleBlockSplittingPolicy(btMultiBodyConstraintSolver* solver);

	/// Destructor
	virtual ~btSingleBlockSplittingPolicy();

	// Documentation inherited
	virtual void split(const btConstraintBlock& blockInput, const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs, btAlignedObjectArray<btConstraintBlock>& blocksOutput);
};

class btDoubleBlockSplittingPolicy : public btBlockSplittingPolicy
{
protected:
	btMultiBodyConstraintSolver* m_solver;

public:
	/// Constructor
	btDoubleBlockSplittingPolicy(btMultiBodyConstraintSolver* solver);

	/// Destructor
	virtual ~btDoubleBlockSplittingPolicy();

	// Documentation inherited
	virtual void split(const btConstraintBlock& blockInput, const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs, btAlignedObjectArray<btConstraintBlock>& blocksOutput);
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
	virtual void solveMultiBodyGroup(btCollisionObject** bodies, int numBodies, btPersistentManifold** manifold, int numManifolds, btTypedConstraint** constraints, int numConstraints, btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer, btDispatcher* dispatcher);

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
	virtual void removeConfig(int configIndex);
};

#endif  // BT_MULTIBODY_BLOCK_CONSTRAINT_SOLVER_H
