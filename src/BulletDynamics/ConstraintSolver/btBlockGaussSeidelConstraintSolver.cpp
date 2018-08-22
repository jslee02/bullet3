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

#include "btBlockGaussSeidelConstraintSolver.h"

#include "BulletDynamics/ConstraintSolver/btConstraintSolver.h"
#include "BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolverMt.h"
#include "BulletDynamics/ConstraintSolver/btSplittingPolicy.h"
#include "LinearMath/btMinMax.h"
#include "LinearMath/btQuickprof.h"

btScalar btBlockGaussSeidelConstraintSolver::solveGroupCacheFriendlySetup(
	btCollisionObject** bodies,
	int numBodies,
	btPersistentManifold** manifoldPtr,
	int numManifolds,
	btTypedConstraint** constraints,
	int numConstraints,
	const btContactSolverInfo& infoGlobal,
	btIDebugDraw* debugDrawer,
	btDispatcher* dispatcher)
{
	m_maxOverrideNumSolverIterations = 0;

	// TODO(JS): Improve
	btSolverConstraintInput input;
	input.m_bodies = bodies;
	input.m_numBodies = numBodies;
	input.m_manifoldPtr = manifoldPtr;
	input.m_numManifolds = numManifolds;
	input.m_constraints = constraints;
	input.m_numConstraints = numConstraints;
	//	block.m_solver = new btSequentialImpulseConstraintSolver();

	if (m_splittingPolicy == 0)
		m_splittingPolicy = new btDummyConstraintBlockGenerator();
	m_blocks.clear();
	bool result = m_splittingPolicy->generate(m_blocks, input);
	btAssert(result);
	// TODO(JS): If result == false, manullay add a simple block with a default
	// constraint solver that takes all the remaining inputs.

	return 0;
}

btScalar btBlockGaussSeidelConstraintSolver::solveGroupCacheFriendlyIterations(
	const btContactSolverInfo& infoGlobal,
	btIDebugDraw* debugDrawer,
	btDispatcher* dispatcher)
{
	BT_PROFILE(
		"btBlockGaussSeidelConstraintSolver::solveGroupCacheFriendlyIterations");

	{
		const int maxIterations =
			m_maxOverrideNumSolverIterations > infoGlobal.m_numIterations
				? m_maxOverrideNumSolverIterations
				: infoGlobal.m_numIterations;

		for (int iteration = 0; iteration < maxIterations; iteration++)
		{
			m_squaredResidual =
				solveSingleIteration(iteration, infoGlobal, debugDrawer, dispatcher);

			if (m_squaredResidual <= infoGlobal.m_leastSquaresResidualThreshold ||
				(iteration >= (maxIterations - 1)))
			{
#ifdef VERBOSE_RESIDUAL_PRINTF
				printf("squared residual = %f at iteration #%d\n", m_squaredResidual,
					   iteration);
#endif
				break;
			}
		}
	}

	return 0;
}

btScalar btBlockGaussSeidelConstraintSolver::solveSingleIteration(
	int iteration, const btContactSolverInfo& infoGlobal,
	btIDebugDraw* debugDrawer, btDispatcher* dispatcher)
{
	btScalar squaredResidual = 0;

	for (int i = 0; i < m_blocks.size(); ++i)
	{
		// Change the sweep direction per every other iteration
		const int index = iteration & 1 ? i : m_blocks.size() - 1 - i;

		const btScalar residual =
			m_blocks[index].solve(infoGlobal, debugDrawer, dispatcher);
		squaredResidual = btMax(squaredResidual, residual * residual);
	}

	return squaredResidual;
}

btScalar btBlockGaussSeidelConstraintSolver::solveGroupCacheFriendlyFinish(
	btCollisionObject** /*bodies*/, int /*numBodies*/,
	const btContactSolverInfo& /*infoGlobal*/)
{
	return 0;
}

btBlockGaussSeidelConstraintSolver::btBlockGaussSeidelConstraintSolver(
	btConstraintBlockGenerator* splittingPolicy)
	: btConstraintSolver(), m_splittingPolicy(splittingPolicy)
{
	// Do nothing
}

btBlockGaussSeidelConstraintSolver::~btBlockGaussSeidelConstraintSolver()
{
	// Do nothing
}

void btBlockGaussSeidelConstraintSolver::setSplittingPolicy(btConstraintBlockGenerator* splittingPolity)
{
	m_splittingPolicy = splittingPolity;
}

btConstraintBlockGenerator* btBlockGaussSeidelConstraintSolver::getSplittingPolicy()
{
	return m_splittingPolicy;
}

const btConstraintBlockGenerator* btBlockGaussSeidelConstraintSolver::getSplittingPolicy() const
{
	return m_splittingPolicy;
}

btScalar btBlockGaussSeidelConstraintSolver::solveGroup(
	btCollisionObject** bodies,
	int numBodies,
	btPersistentManifold** manifold,
	int numManifolds,
	btTypedConstraint** constraints,
	int numConstraints,
	const btContactSolverInfo& infoGlobal,
	btIDebugDraw* debugDrawer,
	btDispatcher* dispatcher)
{
	BT_PROFILE("btBlockGaussSeidelConstraintSolver::solveGroup");

	solveGroupCacheFriendlySetup(bodies, numBodies, manifold, numManifolds,
								 constraints, numConstraints, infoGlobal,
								 debugDrawer, dispatcher);

	solveGroupCacheFriendlyIterations(infoGlobal, debugDrawer, dispatcher);

	solveGroupCacheFriendlyFinish(bodies, numBodies, infoGlobal);

	return 0.f;
}

void btBlockGaussSeidelConstraintSolver::reset()
{
	// Do nothing
}

btConstraintSolverType
btBlockGaussSeidelConstraintSolver::getSolverType() const
{
	return BT_BLOCK_GAUSS_SEIDEL_SOLVER;
}
