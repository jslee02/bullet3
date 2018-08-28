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

#include "btMultiBodyBlockConstraintSolver.h"

#include "LinearMath/btQuickprof.h"

btConstraintBlock::btConstraintBlock()
	: m_constraintConfigId(-1)
{
	// Do nothing
}

btConstraintBlock::btConstraintBlock(
	btTypedConstraint** m_constraints,
	int m_numConstraints,
	btAlignedObjectArray<btSolverBody>* m_solverBodyPool,
	btAlignedObjectArray<btSolverConstraint*>& m_nonContactConstraints,
	btAlignedObjectArray<btSolverConstraint*>& m_normalContactConstraints,
	btAlignedObjectArray<btSolverConstraint*>& m_frictionContactConstraints,
	btAlignedObjectArray<btSolverConstraint*>& m_rollingFrictionContactConstraints,
	btMultiBodyConstraint** m_multiBodyConstraints,
	int m_numMultiBodyConstraints,
	btAlignedObjectArray<btMultiBodySolverConstraint*>& m_multiBodyNonContactConstraints,
	btAlignedObjectArray<btMultiBodySolverConstraint*>& m_multiBodyNormalContactConstraints,
	btAlignedObjectArray<btMultiBodySolverConstraint*>& m_multiBodyFrictionContactConstraints,
	btAlignedObjectArray<btMultiBodySolverConstraint*>& m_multiBodyTorsionalFrictionContactConstraints,
	btMultiBodyJacobianData* m_data)
	: m_constraintConfigId(-1)
{
	// Do nothing
}

btNonSplittingPolicy::btNonSplittingPolicy(btMultiBodyConstraintSolver* solver)
	: m_solver(solver)
{
	// Do nothing
}

btNonSplittingPolicy::~btNonSplittingPolicy()
{
	// Do nothing
}

void btNonSplittingPolicy::split(const btConstraintBlock& blockInput, const btAlignedObjectArray<btBlockConstraintSolverConfig>& availableConfigs, btAlignedObjectArray<btConstraintBlock>& blocksOutput)
{
	btConstraintBlock newBlock = blockInput;
	newBlock.m_solver = m_solver;
	//	newBlock.m_constraints = blockInput.m_constraints;
	//	newBlock.m_numConstraints = blockInput.m_numConstraints;
	//	newBlock.m_multiBodyConstraints = blockInput.m_multiBodyConstraints;
	//	newBlock.m_numMultiBodyConstraints = blockInput.m_numMultiBodyConstraints;

	blocksOutput.push_back(newBlock);
}

btBlockSplittingPolicy::~btBlockSplittingPolicy()
{
	// Do nothing
}

btMultiBodyBlockConstraintSolver::btMultiBodyBlockConstraintSolver()
{
	// Do nothing
}

btMultiBodyBlockConstraintSolver::~btMultiBodyBlockConstraintSolver()
{
	// Do nothing
}

void btMultiBodyBlockConstraintSolver::solveMultiBodyGroup(
	btCollisionObject** bodies,
	int numBodies,
	btPersistentManifold** manifold,
	int numManifolds,
	btTypedConstraint** constraints,
	int numConstraints,
	btMultiBodyConstraint** multiBodyConstraints,
	int numMultiBodyConstraints,
	const btContactSolverInfo& info,
	btIDebugDraw* debugDrawer,
	btDispatcher* dispatcher)
{
	BT_PROFILE("btMultiBodyBlockConstraintSolver::solveMultiBodyGroup");

	// Customized solveMultiBodyGroup() for constraint blocks to avoid calling
	// btSequentialImpulseConstraintSolver::solveGroupConvertConstraints() over and over for each constraint block
	// because the constraint conversion should be done just once by the (parent) block solver.

	solveMultiBodyGroupPrestep(multiBodyConstraints, numMultiBodyConstraints, info, debugDrawer, dispatcher);

	// TODO(JS): Important!! This is just for a test.
	m_multiBodyNonContactConstraints.resize(0);
	m_multiBodyNormalContactConstraints.resize(0);
	m_multiBodyFrictionContactConstraints.resize(0);
	m_multiBodyTorsionalFrictionContactConstraints.resize(0);

	// 1. Setup

	// Convert rigid bodies/multibodies, joints, contacts into constraints.
	solveGroupConvertConstraints(bodies, numBodies, manifold, numManifolds, constraints, numConstraints, info, debugDrawer);

	// TODO(JS): [Debug]
	if (m_tmpSolverBodyPool.size() == 0)
		int a = 10;

	// Split constraints into constraint blocks
	btConstraintBlock input;
	input.m_constraints = constraints;
	input.m_numConstraints = numConstraints;
	input.m_solverBodyPool = &m_tmpSolverBodyPool;
	input.m_nonContactConstraints = m_tmpSolverNonContactConstraintPtrPool;
	input.m_normalContactConstraints = m_tmpSolverContactConstraintPtrPool;
	input.m_frictionContactConstraints = m_tmpSolverContactFrictionConstraintPtrPool;
	input.m_rollingFrictionContactConstraints = m_tmpSolverContactRollingFrictionConstraintPtrPool;
	input.m_multiBodyConstraints = m_tmpMultiBodyConstraints;
	input.m_numMultiBodyConstraints = m_tmpNumMultiBodyConstraints;
	input.m_multiBodyNonContactConstraints = m_multiBodyNonContactConstraintPtrs;
	input.m_multiBodyNormalContactConstraints = m_multiBodyNormalContactConstraintPtrs;
	input.m_multiBodyFrictionContactConstraints = m_multiBodyFrictionContactConstraintPtrs;
	input.m_multiBodyTorsionalFrictionContactConstraints = m_multiBodyTorsionalFrictionContactConstraintPtrs;
	input.m_data = &m_data;
	btAlignedObjectArray<btBlockConstraintSolverConfig> tmp;
	// TODO(JS): This is just for test
	m_splittingPolicy = new btNonSplittingPolicy(new btMultiBodyConstraintSolver());
	btAssert(m_splittingPolicy);
	m_blocks.resize(0);
	m_splittingPolicy->split(input, tmp, m_blocks);

	// TODO(JS): [Debug]
	if (m_tmpSolverBodyPool.size() == 0)
		int a = 10;

	const int maxIterations = m_maxOverrideNumSolverIterations > info.m_numIterations ? m_maxOverrideNumSolverIterations : info.m_numIterations;

	for (int iteration = 0; iteration < maxIterations; ++iteration)
	{
		for (int i = 0; i < m_blocks.size(); ++i)
		{
			// Change the sweep direction per every iteration
			const int index = i & 1 ? i : m_blocks.size() - 1 - i;
			// TODO(JS): Maybe control this by info?

			btConstraintBlock& block = m_blocks[index];
			btMultiBodyConstraintSolver* solver = block.m_solver;
			btAssert(solver);

			// TODO(JS): [Debug]
			if (m_tmpSolverBodyPool.size() == 0)
				int a = 10;

			// Setup constraint blocks
			solver->solveMultiBodyGroupPrestep(
				block.m_multiBodyConstraints, block.m_numMultiBodyConstraints, info, debugDrawer, dispatcher);

			// TODO(JS): [Debug]
			if (m_tmpSolverBodyPool.size() == 0)
				int a = 10;

			// IMPORTANT !!!
			// TODO(JS): Make below function for multibody so that we clear multibody constraint arraies
			// IMPORTANT !!!

			// TODO(JS): [Debug]
			if (m_tmpSolverBodyPool.size() == 0)
				int a = 10;

			solver->solveGroupSolverSpecificInit(
				block.m_solverBodyPool,
				block.m_nonContactConstraints,
				block.m_normalContactConstraints,
				block.m_frictionContactConstraints,
				block.m_rollingFrictionContactConstraints,
				info,
				debugDrawer);

			// TODO(JS): [Debug]
			if (m_tmpSolverBodyPool.size() == 0)
				int a = 10;

			// 2. Iterations: perform Gauss-Seidel iterations
			// TODO(JS): Add split impulse
			m_leastSquaresResidual = solver->solveMultiBodyGroupCacheFriendlyIterationsNew(
				block.m_constraints,
				block.m_numConstraints,
				block.m_solverBodyPool,
				block.m_nonContactConstraints,
				block.m_normalContactConstraints,
				block.m_frictionContactConstraints,
				block.m_rollingFrictionContactConstraints,
				block.m_data,
				block.m_multiBodyNonContactConstraints,
				block.m_multiBodyNormalContactConstraints,
				block.m_multiBodyFrictionContactConstraints,
				block.m_multiBodyTorsionalFrictionContactConstraints,
				info,
				debugDrawer);

			// 3. Finish constraint blocks
			solver->solveMultiBodyGroupCacheFriendlyFinishNew(
				block.m_solverBodyPool,
				block.m_nonContactConstraints,
				block.m_normalContactConstraints,
				block.m_frictionContactConstraints,
				block.m_rollingFrictionContactConstraints,
				block.m_data,
				block.m_multiBodyNonContactConstraints,
				block.m_multiBodyNormalContactConstraints,
				block.m_multiBodyFrictionContactConstraints,
				block.m_multiBodyTorsionalFrictionContactConstraints,
				info);

			solver->solveMultiBodyGroupPoststep(
				block.m_multiBodyConstraints, block.m_numMultiBodyConstraints, info, debugDrawer, dispatcher);
			if (m_leastSquaresResidual <= info.m_leastSquaresResidualThreshold || (iteration >= (maxIterations - 1)))
			{
#ifdef VERBOSE_RESIDUAL_PRINTF
				printf("residual = %f at iteration #%d\n", m_leastSquaresResidual, iteration);
#endif
				break;
			}
		}
	}

	// TODO(JS): For block solver
	// solveMultiBodyGroupCacheFriendlyFinishNew

	solveMultiBodyGroupPoststep(multiBodyConstraints, numMultiBodyConstraints, info, debugDrawer, dispatcher);
}

void btMultiBodyBlockConstraintSolver::setSplittingPolicy(btBlockSplittingPolicy* policy)
{
	m_splittingPolicy = policy;
}

int btMultiBodyBlockConstraintSolver::addConfig(btBlockConstraintSolverConfig& config)
{
	m_configs.push_back(config);
	return m_configs.size();
}

int btMultiBodyBlockConstraintSolver::getNumConfigs() const
{
	return m_configs.size();
}

void btMultiBodyBlockConstraintSolver::removeConfig(int configIndex)
{
	m_configs.removeAtIndex(configIndex);
}
