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

#ifndef BT_CONSTRAINT_SOLVER_FACTORY_H
#define BT_CONSTRAINT_SOLVER_FACTORY_H

#include "BulletDynamics/ConstraintSolver/btConstraintSolver.h"
#include "BulletDynamics/ConstraintSolver/btContactSolverInfo.h"
#include "BulletDynamics/MLCPSolvers/btMLCPSolverInterface.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btThreads.h"

struct btConstraintSolverCreateStructure
{
	btConstraintSolverType m_constraintSolverType;

	/// This is only utilized when \c m_constraintSolverType is
	/// \c BT_MLCP_SOLVER.
	btMLCPSolverType m_mlcpSolverType;

	btConstraintSolverCreateStructure(btConstraintSolverType constraintSolverType, btMLCPSolverType mlcpSolverType = BT_DANTZIG);
};

class btConstraintSolverFactory
{
protected:
public:
	virtual ~btConstraintSolverFactory();

	virtual btConstraintSolver* create(const btConstraintSolverCreateStructure& type) const;
};

#endif  // BT_CONSTRAINT_SOLVER_FACTORY_H
