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

#ifndef BT_SPLITTING_POLICY_H
#define BT_SPLITTING_POLICY_H

#include "BulletDynamics/ConstraintSolver/btConstraintSolver.h"
#include "BulletDynamics/ConstraintSolver/btBlockGaussSeidelConstraintSolver.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btThreads.h"

class btConstraintBlockGenerator
{
public:
	virtual ~btConstraintBlockGenerator();

	/// \return True if all constraint blocks are generated and added to \c block so there is no more input left in
	/// \c input. False, otherise.
	virtual bool generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input) = 0;
};

class btSimpleConstraintBlockGenerator : public btConstraintBlockGenerator
{
protected:
	btConstraintSolver* m_constraintSolver;

public:
	btSimpleConstraintBlockGenerator(btConstraintSolver* constraintSolverType);

	bool generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input);
};

class btContactManifoldConstraintBlockGenerator : public btConstraintBlockGenerator
{
protected:
	btConstraintSolver* m_constraintSolver;

public:
	btContactManifoldConstraintBlockGenerator(btConstraintSolver* constraintSolverType);

	bool generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input);
};

class btDummyConstraintBlockGenerator : public btSimpleConstraintBlockGenerator
{
protected:
public:
	btDummyConstraintBlockGenerator();

	// Documentation inherited
	bool generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input) BT_OVERRIDE;
};

class btChainedBlockGenerators : public btConstraintBlockGenerator
{
protected:
	btAlignedObjectArray<btConstraintBlockGenerator*> m_generatorChain;

public:
	void pushBack(btConstraintBlockGenerator* generator);
	// insert()
	// ...

	bool generate(btAlignedObjectArray<btConstraintBlock>& blocks, btSolverConstraintInput& input) BT_OVERRIDE;
};

#endif  // BT_SPLITTING_POLICY_H
