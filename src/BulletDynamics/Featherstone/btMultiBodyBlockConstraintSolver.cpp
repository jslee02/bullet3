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

btConstraintBlock::btConstraintBlock()
	: m_constraintConfigId(-1)
{
	// Do nothing
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
	// Convert rigid bodies/multibodies, joints, contacts into constraints.
	solveGroupConvertConstraints(bodies, numBodies, manifold, numManifolds, constraints, numConstraints, info, debugDrawer);

	// Split constraints into constraint blocks
	// m_splittingPolicy->split();

	// Setup constraint blocks
	for (int i = 0; i < m_blocks.size(); ++i)
	{
		btMultiBodyConstraintSolver* solver = m_blocks[i].m_solver;

		// Customized solveMultiBodyGroup() for constraint blocks to avoid calling
		// btSequentialImpulseConstraintSolver::solveGroupConvertConstraints() for each constraint blocks because
		// the constraint conversion should be done once by the block solver.
		solver->solveMultiBodyGroupPrestep(
			bodies, numBodies, manifold, numManifolds, constraints, numConstraints, multiBodyConstraints,
			numMultiBodyConstraints, info, debugDrawer, dispatcher);
		solver->solveGroupSolverSpecificInit(info, debugDrawer);
		solver->solveMultiBodyGroupPoststep(
			bodies, numBodies, manifold, numManifolds, constraints, numConstraints, multiBodyConstraints,
			numMultiBodyConstraints, info, debugDrawer, dispatcher);
		// TODO(JS): Consider change the arguments to be similar to blow function calls
	}

	// Perform Gauss-Seidel iterations
	// TODO(JS): Add split impulse
	int maxIterations = m_maxOverrideNumSolverIterations > info.m_numIterations ? m_maxOverrideNumSolverIterations : info.m_numIterations;

	for (int iteration = 0; iteration < maxIterations; ++iteration)
	{
		for (int i = 0; i < m_blocks.size(); ++i)
		{
			// Change the sweep direction per every iteration
			const int index = i & 1 ? i : m_blocks.size() - 1 - i;
			// TODO(JS): Maybe control this by info?

			btConstraintBlock& block = m_blocks[index];
			btMultiBodyConstraintSolver* solver = block.m_solver;
			m_leastSquaresResidual = solver->solveMultiBodySingleIterationNew(
				iteration,
				block.m_constraints,
				block.m_numConstraints,
				block.m_tmpSolverBodyPool,
				block.m_tmpSolverNonContactConstraintPool,
				block.m_tmpSolverContactConstraintPool,
				block.m_tmpSolverContactFrictionConstraintPool,
				block.m_tmpSolverContactRollingFrictionConstraintPool,
				block.m_data,
				block.m_multiBodyNonContactConstraints,
				block.m_multiBodyNormalContactConstraints,
				block.m_multiBodyFrictionContactConstraints,
				block.m_multiBodyTorsionalFrictionContactConstraints,
				info,
				debugDrawer);

			if (m_leastSquaresResidual <= info.m_leastSquaresResidualThreshold || (iteration >= (maxIterations - 1)))
			{
#ifdef VERBOSE_RESIDUAL_PRINTF
				printf("residual = %f at iteration #%d\n", m_leastSquaresResidual, iteration);
#endif
				break;
			}
		}
	}

	// Finialize constraint blocks
	for (int i = 0; i < m_blocks.size(); ++i)
	{
		btConstraintBlock& block = m_blocks[i];
		btMultiBodyConstraintSolver* solver = block.m_solver;

		solver->solveGroupCacheFriendlyFinishNew(
			block.m_tmpSolverBodyPool,
			block.m_tmpSolverNonContactConstraintPool,
			block.m_tmpSolverContactConstraintPool,
			block.m_tmpSolverContactFrictionConstraintPool,
			block.m_tmpSolverContactRollingFrictionConstraintPool,
			block.m_data,
			block.m_multiBodyNonContactConstraints,
			block.m_multiBodyNormalContactConstraints,
			block.m_multiBodyFrictionContactConstraints,
			block.m_multiBodyTorsionalFrictionContactConstraints,
			info);
	}

	// TODO(JS): convert solver body to multibody/rigidbody
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
