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

#include "btSplittingPolicy.h"

#include "LinearMath/btMinMax.h"
#include "LinearMath/btQuickprof.h"
#include "BulletDynamics/ConstraintSolver/btConstraintSolverFactory.h"
#include "BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h"

btConstraintBlock::btConstraintBlock(const btSolverConstraintInput& input, btConstraintSolver* solver)
	: m_input(input), m_solver(solver)
{
	// Do nothing
}

btScalar btConstraintBlock::solve(const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer, btDispatcher* dispatcher)
{
	return m_solver->solveGroup(
		m_input.m_bodies,
		m_input.m_numBodies,
		m_input.m_manifoldPtr,
		m_input.m_numManifolds,
		m_input.m_constraints,
		m_input.m_numConstraints,
		infoGlobal,
		debugDrawer,
		dispatcher);
}

btConstraintBlockGenerator::~btConstraintBlockGenerator()
{
	// Do nothing
}

bool btContactManifoldConstraintBlockGenerator::generate(btAlignedObjectArray<btConstraintBlock> &blocks, btSolverConstraintInput &input)
{
	for (int i = 0; i < input.m_numManifolds; ++i)
	{
		btPersistentManifold* manifold = input.m_manifoldPtr[i];



//		btConstraintBlock block;
//		block.m_solver = m_constraintSolver;
	}
}

btDummyConstraintBlockGenerator::btDummyConstraintBlockGenerator()
	: btSimpleConstraintBlockGenerator(new btSequentialImpulseConstraintSolver())
{
	// Do nothing
}

bool btDummyConstraintBlockGenerator::generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input)
{
	return btSimpleConstraintBlockGenerator::generate(blocks, input);
}

btSimpleConstraintBlockGenerator::btSimpleConstraintBlockGenerator(btConstraintSolver* constraintSolverType)
	: m_constraintSolver(constraintSolverType)
{
	// Do nothing
}

bool btSimpleConstraintBlockGenerator::generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input) BT_OVERRIDE
{
	blocks.clear();
	blocks.push_back(btConstraintBlock(input, m_constraintSolver));

	input.m_numBodies = 0;
	input.m_numManifolds = 0;
	input.m_numConstraints = 0;

	return false;
}

void btChainedBlockGenerators::pushBack(btConstraintBlockGenerator* generator)
{
	m_generatorChain.push_back(generator);
}

bool btChainedBlockGenerators::generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input)
{
	for (int i = 0; i < m_generatorChain.size(); ++i)
	{
		btConstraintBlockGenerator* generator = m_generatorChain[i];
		generator->generate(blocks, input);

		if (input.empty())
			return true;
	}

	return false;
}
