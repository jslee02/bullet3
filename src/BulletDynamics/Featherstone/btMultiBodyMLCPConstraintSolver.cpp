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

void btMultiBodyMLCPConstraintSolver::createMLCPFast(const btContactSolverInfo& infoGlobal)
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
			const btMultiBodySolverConstraint& constraint = *m_allConstraintPtrArray[i];
			const btScalar jacDiag = constraint.m_jacDiagABInv;

			if (!btFuzzyZero(jacDiag))
			{
				// Note that rhsPenetration is currently always zero because the split impulse hasn't been implemented for multibody yet.
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
			if (0)//m_limitDependencies[i]>=0)
			{
				m_lo[i] = -BT_INFINITY;
				m_hi[i] = BT_INFINITY;
			}
			else
			{
				const btMultiBodySolverConstraint& constraint = *m_allConstraintPtrArray[i];
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

			const btMultiBodySolverConstraint& constraint = *m_allConstraintPtrArray[i];

			const btMultiBody* multiBodyA = constraint.m_multiBodyA;
			const btMultiBody* multiBodyB = constraint.m_multiBodyB;

			// TODO(JS): Assumed multibody for now
			assert(multiBodyA);
			assert(multiBodyB);

			btScalar diagA = btScalar(0);

			const btScalar* jacA = &m_data.m_jacobians[constraint.m_jacAindex];
			const btScalar* deltaA = &m_data.m_deltaVelocitiesUnitImpulse[constraint.m_jacAindex];
			const int ndofA  = multiBodyA->getNumDofs() + 6;
			diagA += dot(jacA, deltaA, ndofA);

			const btScalar* jacB = &m_data.m_jacobians[constraint.m_jacBindex];
			const btScalar* deltaB = &m_data.m_deltaVelocitiesUnitImpulse[constraint.m_jacBindex];
			const int ndofB  = multiBodyB->getNumDofs() + 6;
			diagA += dot(jacB, deltaB, ndofB);

			m_A.setElem(i, i, diagA);

			// Computes the off-diagonals of A:
			//   a. The rest of i-th row of A, from A(i, i+1) to A(i, n)
			//   b. The rest of i-th column of A, from A(i+1, i) to A(n, i)
			for (int j = i + 1; j < numConstraints; ++j)
			{
				const btMultiBodySolverConstraint& offDiagConstraint = *m_allConstraintPtrArray[j];

				const btMultiBody* offDiagMultiBodyA = offDiagConstraint.m_multiBodyA;
				const btMultiBody* offDiagMultiBodyB = offDiagConstraint.m_multiBodyB;

				btScalar offDiagA = btScalar(0);

				const btScalar* offDiagJacA = &m_data.m_jacobians[offDiagConstraint.m_jacAindex];
				if (offDiagMultiBodyA == multiBodyA)
				{
					offDiagA += dot(offDiagJacA, deltaA, ndofA);
				}
				else if (offDiagMultiBodyA == multiBodyB)
				{
					offDiagA += dot(offDiagJacA, deltaB, ndofB);
				}

				const btScalar* offDiagJacB = &m_data.m_jacobians[offDiagConstraint.m_jacBindex];
				if (offDiagMultiBodyB == multiBodyA)
				{
					offDiagA += dot(offDiagJacB, deltaA, ndofA);
				}
				else if (offDiagMultiBodyB == multiBodyB)
				{
					offDiagA += dot(offDiagJacB, deltaB, ndofB);
				}

				// Set the off-diagonal values of A. Note that A is symmetric.
				m_A.setElem(i, j, offDiagA);
				m_A.setElem(j, i, offDiagA);
			}
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
				const btMultiBodySolverConstraint& constraint = *m_allConstraintPtrArray[i];
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

bool btMultiBodyMLCPConstraintSolver::solveMLCP(const btContactSolverInfo &infoGlobal)
{
	bool result = true;

	if (m_A.rows()==0)
		return true;

	//if using split impulse, we solve 2 separate (M)LCPs
	if (infoGlobal.m_splitImpulse)
	{
		btMatrixXu Acopy = m_A;
		btAlignedObjectArray<int> limitDependenciesCopy = m_limitDependencies;
//		printf("solve first LCP\n");
		result = m_solver->solveMLCP(m_A, m_b, m_x, m_lo,m_hi, m_limitDependencies,infoGlobal.m_numIterations );
		if (result)
			result = m_solver->solveMLCP(Acopy, m_bSplit, m_xSplit, m_lo,m_hi, limitDependenciesCopy,infoGlobal.m_numIterations );
	}
	else
	{
		result = m_solver->solveMLCP(m_A, m_b, m_x, m_lo,m_hi, m_limitDependencies,infoGlobal.m_numIterations );
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
	btMultiBodyConstraintSolver::solveGroupCacheFriendlySetup(
		bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);

	{
		BT_PROFILE("gather constraint data");

		const int numFrictionPerContact = m_tmpSolverContactConstraintPool.size() == m_tmpSolverContactFrictionConstraintPool.size() ? 1 : 2;

		m_allConstraintPtrArray.resize(0);
		m_limitDependencies.resize(m_multiBodyNonContactConstraints.size()+m_multiBodyNormalContactConstraints.size()+m_multiBodyFrictionContactConstraints.size());
		btAssert(m_limitDependencies.size() == m_multiBodyNonContactConstraints.size()+m_multiBodyNormalContactConstraints.size()+m_multiBodyFrictionContactConstraints.size());

		int dindex = 0;
		for (int i=0;i<m_multiBodyNonContactConstraints.size();++i)
		{
			m_allConstraintPtrArray.push_back(&m_multiBodyNonContactConstraints[i]);
			m_limitDependencies[dindex++] = -1;
		}

		// The btSequentialImpulseConstraintSolver moves all friction constraints at the very end, we can also interleave them instead

		const int firstContactConstraintOffset = dindex;

		if (interleaveContactAndFriction)
		{
			for (int i=0;i<m_multiBodyNormalContactConstraints.size();++i)
			{
				m_allConstraintPtrArray.push_back(&m_multiBodyNormalContactConstraints[i]);
				m_limitDependencies[dindex++] = -1;

				btMultiBodySolverConstraint& frictionContactConstraint1 = m_multiBodyFrictionContactConstraints[i*numFrictionPerContact];

				m_allConstraintPtrArray.push_back(&frictionContactConstraint1);

				const int findex = (frictionContactConstraint1.m_frictionIndex*(1 + numFrictionPerContact));
				m_limitDependencies[dindex++] = findex + firstContactConstraintOffset;

				if (numFrictionPerContact == 2)
				{
					btMultiBodySolverConstraint& frictionContactConstraint2 = m_multiBodyFrictionContactConstraints[i*numFrictionPerContact+1];

					m_allConstraintPtrArray.push_back(&frictionContactConstraint2);
					m_limitDependencies[dindex++] = findex + firstContactConstraintOffset;
				}
			}
		}
		else
		{
			for (int i=0;i<m_multiBodyNormalContactConstraints.size();++i)
			{
				m_allConstraintPtrArray.push_back(&m_multiBodyNormalContactConstraints[i]);
				m_limitDependencies[dindex++] = -1;
			}
			for (int i=0;i<m_multiBodyFrictionContactConstraints.size();++i)
			{
				m_allConstraintPtrArray.push_back(&m_multiBodyFrictionContactConstraints[i]);
				m_limitDependencies[dindex++] = m_multiBodyFrictionContactConstraints[i].m_frictionIndex+firstContactConstraintOffset;
			}
		}

		if (!m_allConstraintPtrArray.size())
		{
			m_A.resize(0,0);
			m_b.resize(0);
			m_x.resize(0);
			m_lo.resize(0);
			m_hi.resize(0);

			return btScalar(0);
		}
	}

	{
		BT_PROFILE("createMLCPFast");
		createMLCPFast(infoGlobal);
	}

	return btScalar(0);
}

btScalar btMultiBodyMLCPConstraintSolver::solveGroupCacheFriendlyIterations(btCollisionObject **bodies, int numBodies, btPersistentManifold **manifoldPtr, int numManifolds, btTypedConstraint **constraints, int numConstraints, const btContactSolverInfo &infoGlobal, btIDebugDraw *debugDrawer)
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
		return btSequentialImpulseConstraintSolver::solveGroupCacheFriendlyIterations(bodies ,numBodies,manifoldPtr, numManifolds,constraints,numConstraints,infoGlobal,debugDrawer);
	}

	// TODO(JS): Add rigidbody--rigidbody and set leastSquaredResidual with the result

	if (result)
	{
		BT_PROFILE("process MLCP results");

		for (int i = 0; i < m_allConstraintPtrArray.size(); ++i)
		{
			btMultiBodySolverConstraint& c = *m_allConstraintPtrArray[i];

			btMultiBody* multiBodyA = c.m_multiBodyA;
			btMultiBody* multiBodyB = c.m_multiBodyB;

			const int ndofA = multiBodyA->getNumDofs() + 6;
			const int ndofB = multiBodyB->getNumDofs() + 6;

			const btScalar deltaImpulse = m_x[i] - c.m_appliedImpulse;
			c.m_appliedPushImpulse = m_x[i];

			applyDeltaVee(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacAindex], deltaImpulse, c.m_deltaVelAindex, ndofA);
#ifdef DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			//note: update of the actual velocities (below) in the multibody does not have to happen now since m_deltaVelocities can be applied after all iterations
			//it would make the multibody solver more like the regular one with m_deltaVelocities being equivalent to btSolverBody::m_deltaLinearVelocity/m_deltaAngularVelocity
			multiBodyA->applyDeltaVeeMultiDof2(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacAindex], deltaImpulse);
#endif //DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			// TODO(JS): RigidBody

			applyDeltaVee(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacBindex], deltaImpulse, c.m_deltaVelBindex, ndofB);
#ifdef DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			//note: update of the actual velocities (below) in the multibody does not have to happen now since m_deltaVelocities can be applied after all iterations
			//it would make the multibody solver more like the regular one with m_deltaVelocities being equivalent to btSolverBody::m_deltaLinearVelocity/m_deltaAngularVelocity
			multiBodyB->applyDeltaVeeMultiDof2(&m_data.m_deltaVelocitiesUnitImpulse[c.m_jacBindex], deltaImpulse);
#endif //DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS
			// TODO(JS): RigidBody

			if (infoGlobal.m_splitImpulse)
			{
				const btScalar deltaImpulse = m_xSplit[i] - c.m_appliedPushImpulse;
				c.m_appliedPushImpulse = m_xSplit[i];

				// TODO(JS): Split impulse
			}
		}
	}

	//return leastSquaredResidual;
	return 0;
}

btMultiBodyMLCPConstraintSolver::btMultiBodyMLCPConstraintSolver(btMLCPSolverInterface *solver)
	: m_solver(solver), m_fallback(0)
{
	// Do nothing
}

btMultiBodyMLCPConstraintSolver::~btMultiBodyMLCPConstraintSolver()
{
	// Do nothing
}

void btMultiBodyMLCPConstraintSolver::setMLCPSolver(btMLCPSolverInterface *solver)
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
