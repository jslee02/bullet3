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

#include "BulletDynamics/Featherstone/btMultiBodyMLCPConstraintSolver.h"

#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraint.h"
#include "BulletDynamics/MLCPSolvers/btMLCPSolverInterface.h"

#define DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS

static bool interleaveContactAndFriction = false;

static btScalar dot(const btScalar* v1, const btScalar* v2, int size)
{
	btScalar result = btScalar(0);
	for (int i = 0; i < size; ++i)
		result += v1[i] * v2[i];

	return result;
}

static btScalar computeRigidBodyConstraintVelocityChange(
	const btVector3& angularDeltaVelocityPerUnitImpulse,
	const btVector3& angularJacobian,
	btScalar invMass)
{
	return angularDeltaVelocityPerUnitImpulse.dot(angularJacobian) + invMass;
}

static btScalar computeRigidBodyConstraintVelocityChange(
	const btVector3& angularDeltaVelocityPerUnitImpulse,
	const btVector3& linearDeltaVelocityPerUnitImpulse,
	const btVector3& linearJacobian,
	const btVector3& angularJacobian)
{
	return angularDeltaVelocityPerUnitImpulse.dot(linearJacobian) + linearDeltaVelocityPerUnitImpulse.dot(angularJacobian);
}

// TODO(JS): Maybe need to chagne the function name to encode this is for diagonal constraints?
static btScalar computeConstraintMatrixElement(
	const btAlignedObjectArray<btSolverBody>& solverBodyPool,
	const btSolverConstraint& constraint)
{
	btScalar ret = btScalar(0);

	const int solverBodyIdA = constraint.m_solverBodyIdA;
	const btSolverBody* solverBodyA = &solverBodyPool[solverBodyIdA];
	ret += computeRigidBodyConstraintVelocityChange(
		constraint.m_relpos1CrossNormal,
		constraint.m_angularComponentA,
		solverBodyA->m_originalBody->getInvMass()
	);

	const int solverBodyIdB = constraint.m_solverBodyIdB;
	const btSolverBody* solverBodyB = &solverBodyPool[solverBodyIdB];
	ret += computeRigidBodyConstraintVelocityChange(
		constraint.m_relpos2CrossNormal,           // angular velocity change
		constraint.m_angularComponentB,            // linear velocity change
		solverBodyB->m_originalBody->getInvMass()  // mass inverse
	);

	return ret;
}

// TODO(JS): Maybe need to chagne the function name to encode this is for off diagonal constraints?
static btScalar computeConstraintMatrixElement(
	const btAlignedObjectArray<btSolverBody>& solverBodyPool,
	const btSolverConstraint& constraint,
	const btSolverConstraint& offDiagConstraint)
{
	btScalar ret = btScalar(0);

	const int solverBodyIdA = constraint.m_solverBodyIdA;
	const int solverBodyIdB = constraint.m_solverBodyIdB;

	const int offDiagSolverBodyIdA = offDiagConstraint.m_solverBodyIdA;
	const int offDiagSolverBodyIdB = offDiagConstraint.m_solverBodyIdB;

	const btSolverBody* offDiagSolverBodyA = &solverBodyPool[offDiagSolverBodyIdA];
	const btSolverBody* offDiagSolverBodyB = &solverBodyPool[offDiagSolverBodyIdB];

	// The constraint associated with offDiagSolverBodyId is bodyA, but different contact point
	if (offDiagSolverBodyIdA == solverBodyIdA)
	{
		// Use the same velocity change per unit impulse of body A, but use the constraint Jacobian of offDiagSolverBodyId
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_angularComponentA,                                                        // angular velocity change
			constraint.m_contactNormal1,                                                           // linear velocity change
			offDiagConstraint.m_relpos1CrossNormal,                                                // angular constraint Jacobian
			offDiagConstraint.m_contactNormal1 * offDiagSolverBodyA->m_originalBody->getInvMass()  // linear constraint Jacobian
		);
	}
	// The constraint associated with offDiagSolverBodyId is bodyB, but different contact point
	else if (offDiagSolverBodyIdA == solverBodyIdB)
	{
		// Use the same velocity change per unit impulse of body B, but use the constraint Jacobian of offDiagSolverBodyId
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_angularComponentB,                                                        // angular velocity change
			constraint.m_contactNormal2,                                                           // linear velocity change
			offDiagConstraint.m_relpos1CrossNormal,                                                // angular constraint Jacobian
			offDiagConstraint.m_contactNormal1 * offDiagSolverBodyA->m_originalBody->getInvMass()  // linear constraint Jacobian
		);
	}

	// The constraint associated with offDiagSolverBodyId is bodyA, but different contact point
	if (offDiagSolverBodyIdB == solverBodyIdA)
	{
		// Use the same velocity change per unit impulse of body A, but use the constraint Jacobian of offDiagSolverBodyId
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_angularComponentA,                                                        // angular velocity change
			constraint.m_contactNormal1,                                                           // linear velocity change
			offDiagConstraint.m_relpos2CrossNormal,                                                // angular constraint Jacobian
			offDiagConstraint.m_contactNormal2 * offDiagSolverBodyB->m_originalBody->getInvMass()  // linear constraint Jacobian
		);
	}
	// The constraint associated with offDiagSolverBodyId is bodyB, but different contact point
	else if (offDiagSolverBodyIdB == solverBodyIdB)
	{
		// Use the same velocity change per unit impulse of body B, but use the constraint Jacobian of offDiagSolverBodyId
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_angularComponentB,                                                        // angular velocity change
			constraint.m_contactNormal2,                                                           // linear velocity change
			offDiagConstraint.m_relpos2CrossNormal,                                                // angular constraint Jacobian
			offDiagConstraint.m_contactNormal2 * offDiagSolverBodyB->m_originalBody->getInvMass()  // linear constraint Jacobian
		);
	}

	return ret;
}

static btScalar computeConstraintMatrixElement(
	const btAlignedObjectArray<btSolverBody>& solverBodyPool,
	const btMultiBodyJacobianData& data,
	const btMultiBodySolverConstraint& constraint)
{
	btScalar ret = btScalar(0);

	const btMultiBody* multiBodyA = constraint.m_multiBodyA;
	const btMultiBody* multiBodyB = constraint.m_multiBodyB;

	if (multiBodyA)
	{
		const btScalar* jacA = &data.m_jacobians[constraint.m_jacAindex];
		const btScalar* deltaA = &data.m_deltaVelocitiesUnitImpulse[constraint.m_jacAindex];
		const int ndofA = multiBodyA->getNumDofs() + 6;
		ret += dot(deltaA, jacA, ndofA);
	}
	else
	{
		const int solverBodyIdA = constraint.m_solverBodyIdA;
		assert(solverBodyIdA != -1);
		const btSolverBody* solverBodyA = &solverBodyPool[solverBodyIdA];
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_relpos1CrossNormal,
			constraint.m_angularComponentA,
			solverBodyA->m_originalBody->getInvMass()
		);
	}

	if (multiBodyB)
	{
		const btScalar* jacB = &data.m_jacobians[constraint.m_jacBindex];
		const btScalar* deltaB = &data.m_deltaVelocitiesUnitImpulse[constraint.m_jacBindex];
		const int ndofB = multiBodyB->getNumDofs() + 6;
		ret += dot(deltaB, jacB, ndofB);
	}
	else
	{
		const int solverBodyIdB = constraint.m_solverBodyIdB;
		assert(solverBodyIdB != -1);
		const btSolverBody* solverBodyB = &solverBodyPool[solverBodyIdB];
		ret += computeRigidBodyConstraintVelocityChange(
			constraint.m_relpos1CrossNormal,
			constraint.m_angularComponentA,
			solverBodyB->m_originalBody->getInvMass()
		);
	}

	return ret;
}

static btScalar computeConstraintMatrixElement(
	const btAlignedObjectArray<btSolverBody>& solverBodyPool,
	const btMultiBodyJacobianData& data,
	const btMultiBodySolverConstraint& constraint,
	const btMultiBodySolverConstraint& offDiagConstraint)
{
	btScalar offDiagA = btScalar(0);

	const btMultiBody* multiBodyA = constraint.m_multiBodyA;
	const btMultiBody* multiBodyB = constraint.m_multiBodyB;

	const btScalar* deltaA = &data.m_deltaVelocitiesUnitImpulse[constraint.m_jacAindex];
	const int ndofA = multiBodyA->getNumDofs() + 6;

	const btScalar* deltaB = &data.m_deltaVelocitiesUnitImpulse[constraint.m_jacBindex];
	const int ndofB = multiBodyB->getNumDofs() + 6;

	const btMultiBody* offDiagMultiBodyA = offDiagConstraint.m_multiBodyA;
	const btMultiBody* offDiagMultiBodyB = offDiagConstraint.m_multiBodyB;

	const btScalar* offDiagJacA = &data.m_jacobians[offDiagConstraint.m_jacAindex];

	if (offDiagMultiBodyA == multiBodyA)
	{
		offDiagA += dot(deltaA, offDiagJacA, ndofA);
	}
	else if (offDiagMultiBodyA == multiBodyB)
	{
		offDiagA += dot(deltaB, offDiagJacA, ndofB);
	}

	const btScalar* offDiagJacB = &data.m_jacobians[offDiagConstraint.m_jacBindex];
	if (offDiagMultiBodyB == multiBodyA)
	{
		offDiagA += dot(deltaA, offDiagJacB, ndofA);
	}
	else if (offDiagMultiBodyB == multiBodyB)
	{
		offDiagA += dot(deltaB, offDiagJacB, ndofB);
	}

	return offDiagA;
}

void btMultiBodyMLCPConstraintSolver::createMLCPFast(const btContactSolverInfo& infoGlobal)
{
	createMLCPFastRigidBody(infoGlobal);
	createMLCPFastMultiBody(infoGlobal);
}

void btMultiBodyMLCPConstraintSolver::createMLCPFastRigidBody(const btContactSolverInfo& infoGlobal)
{
	const int numConstraints = m_allConstraintPtrArray.size();

	// 1. Compute b
	{
		BT_PROFILE("init b (rhs)");

		m_b.resize(numConstraints);
		m_b.setZero();

		m_bSplit.resize(numConstraints);
		m_bSplit.setZero();

		for (int i = 0; i < numConstraints; ++i)
		{
			const btSolverConstraint& constraint = *m_allConstraintPtrArray[i];
			const btScalar jacDiag = constraint.m_jacDiagABInv;

			if (!btFuzzyZero(jacDiag))
			{
				const btScalar rhs = constraint.m_rhs;
				const btScalar rhsPenetration = constraint.m_rhsPenetration;
				m_b[i] = rhs / jacDiag;
				m_bSplit[i] = rhsPenetration / jacDiag;
			}
		}
	}

	// 2. Compute lo and hi
	{
		BT_PROFILE("init lo/ho");

		m_lo.resize(numConstraints);
		m_hi.resize(numConstraints);

		for (int i = 0; i < numConstraints; ++i)
		{
			if (0)  //m_limitDependencies[i]>=0)
			{
				m_lo[i] = -BT_INFINITY;
				m_hi[i] = BT_INFINITY;
			}
			else
			{
				const btSolverConstraint& constraint = *m_allConstraintPtrArray[i];
				m_lo[i] = constraint.m_lowerLimit;
				m_hi[i] = constraint.m_upperLimit;
			}
		}
	}

	// 3. Construct A matrix by using the impulse testing
	{
		BT_PROFILE("Compute A");

		{
			BT_PROFILE("m_A.resize");
			m_A.resize(numConstraints, numConstraints);
		}

		for (int i = 0; i < numConstraints; ++i)
		{
			// Compute the diagonal of A, which is A(i, i)
			const btSolverConstraint& constraint = *m_allConstraintPtrArray[i];
			const btScalar diagA = computeConstraintMatrixElement(m_tmpSolverBodyPool, constraint);
			m_A.setElem(i, i, diagA);

			// TODO(JS): comment is outdated
			// Computes the off-diagonals of A:
			//   a. The rest of i-th row of A, from A(i, i+1) to A(i, n)
			//   b. The rest of i-th column of A, from A(i+1, i) to A(n, i)
			for (int j = i + 1; j < numConstraints; ++j)
			{
				// Set the off-diagonal values of A. Note that A is symmetric.
				const btSolverConstraint& offDiagConstraint = *m_allConstraintPtrArray[j];
				const btScalar offDiagA = computeConstraintMatrixElement(m_tmpSolverBodyPool, constraint, offDiagConstraint);
				m_A.setElem(i, j, offDiagA);
				m_A.setElem(j, i, offDiagA);
			}
		}
	}

	if (1)
	{
		// Add CFM to the diagonal of m_A
		for (int i = 0; i < m_A.rows(); ++i)
		{
			m_A.setElem(i, i, m_A(i, i) + infoGlobal.m_globalCfm / infoGlobal.m_timeStep);
		}
	}

	// 4. Initialize x
	{
		BT_PROFILE("resize/init x");

		m_x.resize(numConstraints);
		m_xSplit.resize(numConstraints);

		if (infoGlobal.m_solverMode & SOLVER_USE_WARMSTARTING)
		{
			for (int i = 0; i < numConstraints; ++i)
			{
				const btSolverConstraint& constraint = *m_allConstraintPtrArray[i];
				m_x[i] = constraint.m_appliedImpulse;
				m_xSplit[i] = constraint.m_appliedPushImpulse;
			}
		}
		else
		{
			m_x.setZero();
			m_xSplit.setZero();
		}
	}
}

void btMultiBodyMLCPConstraintSolver::createMLCPFastMultiBody(const btContactSolverInfo& infoGlobal)
{
	const int multiBodyNumConstraints = m_multiBodyAllConstraintPtrArray.size();

	// 1. Compute b
	{
		BT_PROFILE("init b (rhs)");

		m_multiBodyB.resize(multiBodyNumConstraints);
		m_multiBodyB.setZero();

		m_multiBodyBSplit.resize(multiBodyNumConstraints);
		m_multiBodyBSplit.setZero();

		for (int i = 0; i < multiBodyNumConstraints; ++i)
		{
			const btMultiBodySolverConstraint& constraint = *m_multiBodyAllConstraintPtrArray[i];
			const btScalar jacDiag = constraint.m_jacDiagABInv;

			if (!btFuzzyZero(jacDiag))
			{
				// Note that rhsPenetration is currently always zero because the split impulse hasn't been implemented for multibody yet.
				const btScalar rhs = constraint.m_rhs;
				const btScalar rhsPenetration = constraint.m_rhsPenetration;
				m_multiBodyB[i] = rhs / jacDiag;
				m_multiBodyBSplit[i] = rhsPenetration / jacDiag;
			}
		}
	}

	// 2. Compute lo and hi
	{
		BT_PROFILE("init lo/ho");

		m_multiBodyLo.resize(multiBodyNumConstraints);
		m_multiBodyHi.resize(multiBodyNumConstraints);

		for (int i = 0; i < multiBodyNumConstraints; ++i)
		{
			if (0)  //m_limitDependencies[i]>=0)
			{
				m_multiBodyLo[i] = -BT_INFINITY;
				m_multiBodyHi[i] = BT_INFINITY;
			}
			else
			{
				const btMultiBodySolverConstraint& constraint = *m_multiBodyAllConstraintPtrArray[i];
				m_multiBodyLo[i] = constraint.m_lowerLimit;
				m_multiBodyHi[i] = constraint.m_upperLimit;
			}
		}
	}

	// 3. Construct A matrix by using the impulse testing
	{
		BT_PROFILE("Compute A");

		{
			BT_PROFILE("m_A.resize");
			m_multiBodyA.resize(multiBodyNumConstraints, multiBodyNumConstraints);
		}

		for (int i = 0; i < multiBodyNumConstraints; ++i)
		{
			// Compute the diagonal of A, which is A(i, i)
			const btMultiBodySolverConstraint& constraint = *m_multiBodyAllConstraintPtrArray[i];
			const btScalar diagA = computeConstraintMatrixElement(m_tmpSolverBodyPool, m_data, constraint);
			m_multiBodyA.setElem(i, i, diagA);

			// Computes the off-diagonals of A:
			//   a. The rest of i-th row of A, from A(i, i+1) to A(i, n)
			//   b. The rest of i-th column of A, from A(i+1, i) to A(n, i)
			for (int j = i + 1; j < multiBodyNumConstraints; ++j)
			{
				const btMultiBodySolverConstraint& offDiagConstraint = *m_multiBodyAllConstraintPtrArray[j];
				const btScalar offDiagA = computeConstraintMatrixElement(m_tmpSolverBodyPool, m_data, constraint, offDiagConstraint);

				// Set the off-diagonal values of A. Note that A is symmetric.
				m_multiBodyA.setElem(i, j, offDiagA);
				m_multiBodyA.setElem(j, i, offDiagA);
			}
		}
	}

	if (1)
	{
		// Add CFM to the diagonal of m_A
		for (int i = 0; i < m_multiBodyA.rows(); ++i)
		{
			m_multiBodyA.setElem(i, i, m_multiBodyA(i, i) + infoGlobal.m_globalCfm / infoGlobal.m_timeStep);
		}
	}

	// 4. Initialize x
	{
		BT_PROFILE("resize/init x");

		m_multiBodyX.resize(multiBodyNumConstraints);
		m_multiBodyXSplit.resize(multiBodyNumConstraints);

		if (infoGlobal.m_solverMode & SOLVER_USE_WARMSTARTING)
		{
			for (int i = 0; i < multiBodyNumConstraints; ++i)
			{
				const btMultiBodySolverConstraint& constraint = *m_multiBodyAllConstraintPtrArray[i];
				m_multiBodyX[i] = constraint.m_appliedImpulse;
				m_multiBodyXSplit[i] = constraint.m_appliedPushImpulse;
			}
		}
		else
		{
			m_multiBodyX.setZero();
			m_multiBodyXSplit.setZero();
		}
	}
}

bool btMultiBodyMLCPConstraintSolver::solveMLCP(const btContactSolverInfo& infoGlobal)
{
	bool result = true;

	if (m_A.rows() != 0)
	{
		// If using split impulse, we solve 2 separate (M)LCPs
		if (infoGlobal.m_splitImpulse)
		{
			const btMatrixXu Acopy = m_multiBodyA;
			const btAlignedObjectArray<int> limitDependenciesCopy = m_limitDependencies;
			// TODO(JS): Do we really need these copies when solveMLCP takes them as const?

			result = m_solver->solveMLCP(m_A, m_b, m_x, m_lo, m_hi, m_limitDependencies, infoGlobal.m_numIterations);
			if (result)
				result = m_solver->solveMLCP(Acopy, m_bSplit, m_xSplit, m_lo, m_hi, limitDependenciesCopy, infoGlobal.m_numIterations);
		}
		else
		{
			result = m_solver->solveMLCP(m_A, m_b, m_x, m_lo, m_hi, m_limitDependencies, infoGlobal.m_numIterations);
		}
	}

	if (!result)
		return false;

	if (m_multiBodyA.rows() != 0)
	{
		// If using split impulse, we solve 2 separate (M)LCPs
		if (infoGlobal.m_splitImpulse)
		{
			const btMatrixXu Acopy = m_multiBodyA;
			const btAlignedObjectArray<int> limitDependenciesCopy = m_multiBodyLimitDependencies;
			// TODO(JS): Do we really need these copies when solveMLCP takes them as const?

			result = m_solver->solveMLCP(m_multiBodyA, m_multiBodyB, m_multiBodyX, m_multiBodyLo, m_multiBodyHi, m_multiBodyLimitDependencies, infoGlobal.m_numIterations);
			if (result)
				result = m_solver->solveMLCP(Acopy, m_multiBodyBSplit, m_multiBodyXSplit, m_multiBodyLo, m_multiBodyHi, limitDependenciesCopy, infoGlobal.m_numIterations);
		}
		else
		{
			result = m_solver->solveMLCP(m_multiBodyA, m_multiBodyB, m_multiBodyX, m_multiBodyLo, m_multiBodyHi, m_multiBodyLimitDependencies, infoGlobal.m_numIterations);
		}
	}

	return result;
}

btScalar btMultiBodyMLCPConstraintSolver::solveGroupCacheFriendlySetup(
	btCollisionObject** bodies,
	int numBodies,
	btPersistentManifold** manifoldPtr,
	int numManifolds,
	btTypedConstraint** constraints,
	int numConstraints,
	const btContactSolverInfo& infoGlobal,
	btIDebugDraw* debugDrawer)
{
	printf("MLCP\n");

	// 1. Setup for rigid-bodies
	btMultiBodyConstraintSolver::solveGroupCacheFriendlySetup(
		bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);

	// 2. Setup for multi-bodies
	//   a. Collect all different kinds of constraint as pointers into one array, m_allConstraintPtrArray
	//   b. Set the index array for frictional contact constraints, m_limitDependencies
	{
		BT_PROFILE("gather constraint data");

		const int numFrictionPerContact = m_tmpSolverContactConstraintPool.size() == m_tmpSolverContactFrictionConstraintPool.size() ? 1 : 2;
		const int numtiBodyNumFrictionPerContact = m_multiBodyNormalContactConstraints.size() == m_multiBodyFrictionContactConstraints.size() ? 1 : 2;

		m_allConstraintPtrArray.resize(0);
		m_multiBodyAllConstraintPtrArray.resize(0);

		m_multiBodyLimitDependencies.resize(m_tmpSolverNonContactConstraintPool.size() + m_tmpSolverContactConstraintPool.size() + m_tmpSolverContactFrictionConstraintPool.size() + m_multiBodyNonContactConstraints.size() + m_multiBodyNormalContactConstraints.size() + m_multiBodyFrictionContactConstraints.size());
		btAssert(m_multiBodyLimitDependencies.size() == m_tmpSolverNonContactConstraintPool.size() + m_tmpSolverContactConstraintPool.size() + m_tmpSolverContactFrictionConstraintPool.size() + m_multiBodyNonContactConstraints.size() + m_multiBodyNormalContactConstraints.size() + m_multiBodyFrictionContactConstraints.size());

		int dindex = 0;

		// i. Setup for rigid bodies

		for (int i = 0; i < m_tmpSolverNonContactConstraintPool.size(); ++i)
		{
			m_allConstraintPtrArray.push_back(&m_tmpSolverNonContactConstraintPool[i]);
			m_multiBodyLimitDependencies[dindex++] = -1;
		}

		// The btSequentialImpulseConstraintSolver moves all friction constraints at the very end, we can also interleave them instead

		int firstContactConstraintOffset = dindex;

		if (interleaveContactAndFriction)
		{
			for (int i = 0; i < m_tmpSolverContactConstraintPool.size(); i++)
			{
				m_allConstraintPtrArray.push_back(&m_tmpSolverContactConstraintPool[i]);
				m_multiBodyLimitDependencies[dindex++] = -1;
				m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i * numFrictionPerContact]);
				int findex = (m_tmpSolverContactFrictionConstraintPool[i * numFrictionPerContact].m_frictionIndex * (1 + numFrictionPerContact));
				m_multiBodyLimitDependencies[dindex++] = findex + firstContactConstraintOffset;
				if (numFrictionPerContact == 2)
				{
					m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i * numFrictionPerContact + 1]);
					m_multiBodyLimitDependencies[dindex++] = findex + firstContactConstraintOffset;
				}
			}
		}
		else
		{
			for (int i = 0; i < m_tmpSolverContactConstraintPool.size(); i++)
			{
				m_allConstraintPtrArray.push_back(&m_tmpSolverContactConstraintPool[i]);
				m_multiBodyLimitDependencies[dindex++] = -1;
			}
			for (int i = 0; i < m_tmpSolverContactFrictionConstraintPool.size(); i++)
			{
				m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i]);
				m_multiBodyLimitDependencies[dindex++] = m_tmpSolverContactFrictionConstraintPool[i].m_frictionIndex + firstContactConstraintOffset;
			}
		}

		// ii. Setup for multibodies

		for (int i = 0; i < m_multiBodyNonContactConstraints.size(); ++i)
		{
			m_multiBodyAllConstraintPtrArray.push_back(&m_multiBodyNonContactConstraints[i]);
			m_multiBodyLimitDependencies[dindex++] = -1;
		}

		// The btSequentialImpulseConstraintSolver moves all friction constraints at the very end, we can also interleave them instead

		firstContactConstraintOffset = dindex;

		if (interleaveContactAndFriction)
		{
			for (int i = 0; i < m_multiBodyNormalContactConstraints.size(); ++i)
			{
				m_multiBodyAllConstraintPtrArray.push_back(&m_multiBodyNormalContactConstraints[i]);
				m_multiBodyLimitDependencies[dindex++] = -1;

				btMultiBodySolverConstraint& frictionContactConstraint1 = m_multiBodyFrictionContactConstraints[i * numtiBodyNumFrictionPerContact];
				m_multiBodyAllConstraintPtrArray.push_back(&frictionContactConstraint1);

				const int findex = (frictionContactConstraint1.m_frictionIndex * (1 + numtiBodyNumFrictionPerContact)) + firstContactConstraintOffset;

				m_multiBodyLimitDependencies[dindex++] = findex;

				if (numtiBodyNumFrictionPerContact == 2)
				{
					btMultiBodySolverConstraint& frictionContactConstraint2 = m_multiBodyFrictionContactConstraints[i * numtiBodyNumFrictionPerContact + 1];
					m_multiBodyAllConstraintPtrArray.push_back(&frictionContactConstraint2);

					m_multiBodyLimitDependencies[dindex++] = findex;
				}
			}
		}
		else
		{
			for (int i = 0; i < m_multiBodyNormalContactConstraints.size(); ++i)
			{
				m_multiBodyAllConstraintPtrArray.push_back(&m_multiBodyNormalContactConstraints[i]);
				m_multiBodyLimitDependencies[dindex++] = -1;
			}
			for (int i = 0; i < m_multiBodyFrictionContactConstraints.size(); ++i)
			{
				m_multiBodyAllConstraintPtrArray.push_back(&m_multiBodyFrictionContactConstraints[i]);
				m_multiBodyLimitDependencies[dindex++] = m_multiBodyFrictionContactConstraints[i].m_frictionIndex + firstContactConstraintOffset;
			}
		}

		if (!m_multiBodyAllConstraintPtrArray.size())
		{
			m_multiBodyA.resize(0, 0);
			m_multiBodyB.resize(0);
			m_multiBodyX.resize(0);
			m_multiBodyLo.resize(0);
			m_multiBodyHi.resize(0);

			return btScalar(0);
		}
	}

	// Construct MLCP terms
	{
		BT_PROFILE("createMLCPFast");
		createMLCPFast(infoGlobal);
	}

	return btScalar(0);
}

btScalar btMultiBodyMLCPConstraintSolver::solveGroupCacheFriendlyIterations(btCollisionObject** bodies, int numBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer)
{
	bool result = true;
	{
		BT_PROFILE("solveMLCP");
		result = solveMLCP(infoGlobal);
	}

	// Fallback to btSequentialImpulseConstraintSolver::solveGroupCacheFriendlyIterations if the solution isn't valid.
	if (!result)
	{
		m_fallback++;
		return btMultiBodyConstraintSolver::solveGroupCacheFriendlyIterations(bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);
	}

	// TODO(JS): Add rigidbody--rigidbody and set leastSquaredResidual with the result

	if (result)
	{
		BT_PROFILE("process MLCP results");

		for (int i = 0; i < m_allConstraintPtrArray.size(); ++i)
		{
			// TODO(JS): Not implemented.
		}

		for (int i = 0; i < m_multiBodyAllConstraintPtrArray.size(); ++i)
		{
			btMultiBodySolverConstraint& c = *m_multiBodyAllConstraintPtrArray[i];

			btMultiBody* multiBodyA = c.m_multiBodyA;
			btMultiBody* multiBodyB = c.m_multiBodyB;

			const int ndofA = multiBodyA->getNumDofs() + 6;
			const int ndofB = multiBodyB->getNumDofs() + 6;

			const btScalar deltaImpulse = m_multiBodyX[i] - c.m_appliedImpulse;
			c.m_appliedPushImpulse = m_multiBodyX[i];

			applyDeltaVee(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacAindex], deltaImpulse, c.m_deltaVelAindex, ndofA);
#ifdef DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			//note: update of the actual velocities (below) in the multibody does not have to happen now since m_deltaVelocities can be applied after all iterations
			//it would make the multibody solver more like the regular one with m_deltaVelocities being equivalent to btSolverBody::m_deltaLinearVelocity/m_deltaAngularVelocity
			multiBodyA->applyDeltaVeeMultiDof2(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacAindex], deltaImpulse);
#endif  //DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS \
	// TODO(JS): RigidBody

			applyDeltaVee(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacBindex], deltaImpulse, c.m_deltaVelBindex, ndofB);
#ifdef DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			//note: update of the actual velocities (below) in the multibody does not have to happen now since m_deltaVelocities can be applied after all iterations
			//it would make the multibody solver more like the regular one with m_deltaVelocities being equivalent to btSolverBody::m_deltaLinearVelocity/m_deltaAngularVelocity
			multiBodyB->applyDeltaVeeMultiDof2(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacBindex], deltaImpulse);
#endif  //DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS \
	// TODO(JS): RigidBody

			if (infoGlobal.m_splitImpulse)
			{
				const btScalar deltaImpulse = m_multiBodyXSplit[i] - c.m_appliedPushImpulse;
				c.m_appliedPushImpulse = m_multiBodyXSplit[i];

				// TODO(JS): Split impulse
			}
		}
	}

	//return leastSquaredResidual;
	return 0;
}

btMultiBodyMLCPConstraintSolver::btMultiBodyMLCPConstraintSolver(btMLCPSolverInterface* solver)
	: m_solver(solver), m_fallback(0)
{
	// Do nothing
}

btMultiBodyMLCPConstraintSolver::~btMultiBodyMLCPConstraintSolver()
{
	// Do nothing
}

void btMultiBodyMLCPConstraintSolver::setMLCPSolver(btMLCPSolverInterface* solver)
{
	m_solver = solver;
}

int btMultiBodyMLCPConstraintSolver::getNumFallbacks() const
{
	return m_fallback;
}

void btMultiBodyMLCPConstraintSolver::setNumFallbacks(int num)
{
	m_fallback = num;
}

btConstraintSolverType btMultiBodyMLCPConstraintSolver::getSolverType() const
{
	return BT_MLCP_SOLVER;
}
