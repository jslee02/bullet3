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

#include "btConstraintSolverFactory.h"

#include "BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h"
#include "BulletDynamics/MLCPSolvers/btMLCPSolver.h"
#include "BulletDynamics/MLCPSolvers/btMLCPSolverFactory.h"
#include "LinearMath/btMinMax.h"
#include "LinearMath/btQuickprof.h"

btConstraintSolverCreateStructure::btConstraintSolverCreateStructure(btConstraintSolverType constraintSolverType, btMLCPSolverType mlcpSolverType)
	: m_constraintSolverType(constraintSolverType), m_mlcpSolverType(mlcpSolverType)
{
	// Do nothing
}

btConstraintSolverFactory::~btConstraintSolverFactory()
{
	// Do nothing
}

btConstraintSolver* btConstraintSolverFactory::create(const btConstraintSolverCreateStructure& input) const
{
	if (input.m_constraintSolverType == BT_SEQUENTIAL_IMPULSE_SOLVER)
	{
		return new btSequentialImpulseConstraintSolver();
	}
	else if (input.m_constraintSolverType == BT_MLCP_SOLVER)
	{
		btMLCPSolverFactory factory;
		btMLCPSolverInterface* mlcpSolver = factory.create(input.m_mlcpSolverType);
		return new btMLCPSolver(mlcpSolver);
	}
	else
	{
		// TODO(JS): Populate
		return 0;
	}
}
