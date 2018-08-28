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
	btDispatcher* /*dispatcher*/)
{
	m_tmpMultiBodyConstraints = multiBodyConstraints;
	m_tmpNumMultiBodyConstraints = numMultiBodyConstraints;

	// 1. Convert rigid bodies/multibodies, joints, contacts into constraints.
	solveGroupCacheFriendlySetup(bodies, numBodies, manifold, numManifolds, constraints, numConstraints, info, debugDrawer);

	// 2. Split constraints into constraint blocks
	btConstraintBlock input;
	input.m_constraints = constraints;
	input.m_numConstraints = numConstraints;
	input.m_solverBodyPool = &m_tmpSolverBodyPool;
	input.m_nonContactConstraints = m_tmpSolverNonContactConstraintPool;
	input.m_normalContactConstraints = m_tmpSolverContactConstraintPool;
	input.m_frictionContactConstraints = m_tmpSolverContactFrictionConstraintPool;
	input.m_rollingFrictionContactConstraints = m_tmpSolverContactRollingFrictionConstraintPool;
	input.m_multiBodyConstraints = m_tmpMultiBodyConstraints;
	input.m_numMultiBodyConstraints = m_tmpNumMultiBodyConstraints;
	input.m_multiBodyNonContactConstraints = m_multiBodyNonContactConstraints;
	input.m_multiBodyNormalContactConstraints = m_multiBodyNormalContactConstraints;
	input.m_multiBodyFrictionContactConstraints = m_multiBodyFrictionContactConstraints;
	input.m_multiBodyTorsionalFrictionContactConstraints = m_multiBodyTorsionalFrictionContactConstraints;
	input.m_data = &m_data;

	input.m_bodies = bodies;
	input.m_numBodies = numBodies;
	input.m_manifold = manifold;
	input.m_numManifolds = numManifolds;

	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_constraints = constraints;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_numConstraints = numConstraints;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_solverBodyPool = &m_tmpSolverBodyPool;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_nonContactConstraints = m_tmpSolverNonContactConstraintPool;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_normalContactConstraints = m_tmpSolverContactConstraintPool;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_frictionContactConstraints = m_tmpSolverContactFrictionConstraintPool;
	input.m_multiBodyConstraintSet.m_rigidBodyConstraints.m_rollingFrictionContactConstraints = m_tmpSolverContactRollingFrictionConstraintPool;

	input.m_multiBodyConstraintSet.m_multiBodyConstraints = multiBodyConstraints;
	input.m_multiBodyConstraintSet.m_numMultiBodyConstraints = numMultiBodyConstraints;
	input.m_multiBodyConstraintSet.m_multiBodyNonContactConstraints = m_multiBodyNonContactConstraints;
	input.m_multiBodyConstraintSet.m_multiBodyNormalContactConstraints = m_multiBodyNormalContactConstraints;
	input.m_multiBodyConstraintSet.m_multiBodyFrictionContactConstraints = m_multiBodyFrictionContactConstraints;
	input.m_multiBodyConstraintSet.m_multiBodyTorsionalFrictionContactConstraints = m_multiBodyTorsionalFrictionContactConstraints;
	input.m_multiBodyConstraintSet.m_data = &m_data;

	btAlignedObjectArray<btBlockConstraintSolverConfig> tmp;
	// TODO(JS): This is just for test
	m_splittingPolicy = new btNonSplittingPolicy(new btMultiBodyConstraintSolver());
	btAssert(m_splittingPolicy);
	m_blocks.resize(0);
	m_splittingPolicy->split(input, tmp, m_blocks);

	// 3. Setup constraint solvers
	for (int i = 0; i < m_blocks.size(); ++i)
	{
		btConstraintBlock& block = m_blocks[i];
		btMultiBodyConstraintSolver* solver = block.m_solver;
		btAssert(solver);
		solver->solveGroupConvertConstraintPrestep(block.m_bodies, block.m_numBodies, block.m_manifold, block.m_numManifolds, constraints, block.m_numConstraints, info, debugDrawer);
		solver->setMultiBodyConstraints(block.m_multiBodyConstraintSet);
		solver->solveGroupConvertConstraintPoststep(block.m_bodies, block.m_numBodies, block.m_manifold, block.m_numManifolds, constraints, block.m_numConstraints, info, debugDrawer);

		// refine constraint blocks

		// assign constraint blocks

		// write back constraint blocks
	}

	// 4. Gauss-Seidel iterations

	const int maxIterations = m_maxOverrideNumSolverIterations > info.m_numIterations ? m_maxOverrideNumSolverIterations : info.m_numIterations;

	m_leastSquaresResidual = 0;

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

			// Setup constraint blocks

			// 2. Iterations: perform Gauss-Seidel iterations
			// TODO(JS): Add split impulse
			btScalar newSquaredResidual = solver->solveGroupCacheFriendlyIterations(block.m_bodies, block.m_numBodies, block.m_manifold, block.m_numManifolds, block.m_constraints, block.m_numConstraints, info, debugDrawer);

			m_leastSquaresResidual = btMax(m_leastSquaresResidual, newSquaredResidual);
		}

		if (m_leastSquaresResidual <= info.m_leastSquaresResidualThreshold || (iteration >= (maxIterations - 1)))
		{
#ifdef VERBOSE_RESIDUAL_PRINTF
			printf("residual = %f at iteration #%d\n", m_leastSquaresResidual, iteration);
#endif
			break;
		}
	}

	// 5. Finish
	for (int i = 0; i < m_blocks.size(); ++i)
	{
		btConstraintBlock& block = m_blocks[i];
		btMultiBodyConstraintSolver* solver = block.m_solver;
		btAssert(solver);

		// 3. Finish constraint blocks
		m_leastSquaresResidual = solver->solveGroupCacheFriendlyFinish(block.m_bodies, block.m_numBodies, info);
	}

	m_tmpMultiBodyConstraints = 0;
	m_tmpNumMultiBodyConstraints = 0;
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
