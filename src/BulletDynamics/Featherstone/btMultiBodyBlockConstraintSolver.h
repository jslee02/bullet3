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

struct btBlockConstraintSolverConfig
{
	int m_solverType;  //SI or MLCP Dantzig
	//to be decided: full or subset of

	btContactSolverInfo m_info;
};

struct btConstraintBlock
{
	btConstraintBlock()
		: m_constraintConfigId(-1)
	{
		// Do nothing
	}

	btCollisionObject** m_bodies;
	int m_numBodies;
	btPersistentManifold** m_manifold;
	int m_numManifolds;
	btTypedConstraint** m_constraints;
	int m_numConstraints;
	btMultiBodyConstraint** m_multiBodyConstraints;
	int m_numMultiBodyConstraints;

	int m_constraintConfigId;
};

class btBlockSplittingPolicy
{
protected:
	virtual ~btBlockSplittingPolicy()
	{
		// Do nothing
	}

public:
	virtual void split(
		const btConstraintBlock& blockInput,
		const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs,
		btAlignedObjectArray<btConstraintBlock>& blocksOutput) = 0;
};

class btMultiBodyBlockConstraintSolver : public btMultiBodyConstraintSolver
{
public:
	btMultiBodyBlockConstraintSolver();
	virtual ~btMultiBodyBlockConstraintSolver();

	virtual void solveMultiBodyGroup(btCollisionObject** bodies, int numBodies, btPersistentManifold** manifold, int numManifolds, btTypedConstraint** constraints, int numConstraints, btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer, btDispatcher* dispatcher);

	virtual void setSplittingPolicy(btBlockSplittingPolicy* policy);

	virtual int addConfig(btBlockConstraintSolverConfig& config);
	virtual int getNumConfigs() const;
	virtual void removeConfig(int configIndex);  //in range 0..numConfigs
};

#endif  //BT_MULTIBODY_BLOCK_CONSTRAINT_SOLVER_H
