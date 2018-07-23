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

static bool interleaveContactAndFriction = false;

btScalar dot(const btScalar* v1, const btScalar* v2, int size)
{
	btScalar result = btScalar(0);
	for (int i = 0; i < size; ++i)
		result += v1[i] * v2[i];

	return result;
}

void print(const btMatrixXu& mat)
{
	const int rows = mat.rows();
	const int cols = mat.cols();

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			printf("%5.2f", mat(i, j));

			if (j != cols - 1)
				printf(", ");
		}
		printf("\n");
	}
}

// TODO(JS): Consider more constraint types other than normal contact constraints
// TODO(JS): Handle btRigidBody cases
// TODO(JS): Consider self collision case
void btMultiBodyMLCPConstraintSolver::createMLCPFast(const btContactSolverInfo& infoGlobal)
{
	printf("btMultiBodyMLCPConstraintSolver::createMLCPFast\n");

	// TODO(JS): Generalize this other than normal contact constraints
	const btMultiBodyConstraintArray& allConstraints = m_multiBodyNormalContactConstraints;

	const int numContactRows = interleaveContactAndFriction ? 3 : 1;

	// 1. Compute b
	// TODO(JS): Change normal contact constraints to m_allConstraintPtrArray
	const int numConstraintRows = allConstraints.size();
	const int n = numConstraintRows;
	{
		BT_PROFILE("init b (rhs)");

		m_b.resize(numConstraintRows);
		m_b.setZero();

		m_bSplit.resize(numConstraintRows);
		m_bSplit.setZero();

		for (int i = 0; i < numConstraintRows; ++i)
		{
			const btMultiBodySolverConstraint& constraint = allConstraints[i];
			const btScalar jacDiag = constraint.m_jacDiagABInv;
			if (!btFuzzyZero(jacDiag))
			{
				const btScalar rhs = constraint.m_rhs;
				const btScalar rhsPenetration = constraint.m_rhsPenetration;
				m_b[i] = rhs / jacDiag;
				m_bSplit[i] = rhsPenetration / jacDiag;

				// Note that rhsPenetration is currently always zero because the split impulse hasn't been implemented for multibody yet.
			}
		}
	}

	// 2. Compute lo and hi
	m_lo.resize(numConstraintRows);
	m_hi.resize(numConstraintRows);
	{
		BT_PROFILE("init lo/ho");

		for (int i = 0; i < numConstraintRows; ++i)
		{
			if (0)//m_limitDependencies[i]>=0)
			{
				m_lo[i] = -BT_INFINITY;
				m_hi[i] = BT_INFINITY;
			}
			else
			{
				const btMultiBodySolverConstraint& constraint = allConstraints[i];
				m_lo[i] = constraint.m_lowerLimit;
				m_hi[i] = constraint.m_upperLimit;
			}
		}
	}

	// TODO(JS): Maybe we could merge the above two for loops as one.

	// 3. Construct A matrix by using the impulse testing
	{
		BT_PROFILE("Compute A");

		// TODO(JS): Use m_allConstraintPtrArray
		// const int numConstraintRows = m_allConstraintPtrArray.size();
		const int numConstraints = allConstraints.size();

		{
			BT_PROFILE("m_A.resize");
			m_A.resize(numConstraints, numConstraints);
		}

		// Construct A only for normal contact constraints
		for (int i = 0; i < allConstraints.size(); ++i)
		{
			// Constraint for the diagonal of A where the index is (i, i)
			const btMultiBodySolverConstraint& constraint = allConstraints[i];


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

			for (int j = i + 1; j < allConstraints.size(); ++j)
			{
				// Constraint for the off-diagonal of A where the index is (i, j)
				const btMultiBodySolverConstraint& offDiagConstraint = m_multiBodyNormalContactConstraints[j];

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

				// Set off-diagonal elements, which are of symmetric matrix A
				m_A.setElem(i, j, offDiagA);
				m_A.setElem(j, i, offDiagA);
			}
		}
	}

	// 4. Initialize x
	{
		BT_PROFILE("resize/init x");

		m_x.resize(numConstraintRows);
		m_xSplit.resize(numConstraintRows);

		if (infoGlobal.m_solverMode & SOLVER_USE_WARMSTARTING)
		{
			for (int i = 0; i < m_allConstraintPtrArray.size(); ++i)
			{
				const btMultiBodySolverConstraint& constraint = allConstraints[i];
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

	// TODO(JS): Just for quick test. Must be removed.
	m_limitDependencies.resize(numConstraintRows);
	for (int i = 0; i < numConstraintRows; ++i)
		m_limitDependencies[i] = -1;
}

bool btMultiBodyMLCPConstraintSolver::solveMLCP(const btContactSolverInfo &infoGlobal)
{
	bool success = true;

	if (m_A.rows()==0)
		return true;

	//if using split impulse, we solve 2 separate (M)LCPs
	if (infoGlobal.m_splitImpulse)
	{
		btMatrixXu Acopy = m_A;
		btAlignedObjectArray<int> limitDependenciesCopy = m_limitDependencies;
//		printf("solve first LCP\n");
		success = m_solver->solveMLCP(m_A, m_b, m_x, m_lo,m_hi, m_limitDependencies,infoGlobal.m_numIterations );
		if (success)
			success = m_solver->solveMLCP(Acopy, m_bSplit, m_xSplit, m_lo,m_hi, limitDependenciesCopy,infoGlobal.m_numIterations );
	}
	else
	{
		success = m_solver->solveMLCP(m_A, m_b, m_x, m_lo,m_hi, m_limitDependencies,infoGlobal.m_numIterations );
	}

	return success;
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
//	printf("btMultiBodyMLCPConstraintSolver::solveGroupCacheFriendlySetup\n");

	m_multiBodyNormalContactConstraints.resize(0);

	m_data.m_jacobians.resize(0);
	m_data.m_deltaVelocitiesUnitImpulse.resize(0);
	m_data.m_deltaVelocities.resize(0);

	for (int i = 0; i < numBodies; ++i)
	{
		const btMultiBodyLinkCollider* fc = btMultiBodyLinkCollider::upcast(bodies[i]);
		if (fc)
		{
			fc->m_multiBody->setCompanionId(-1);
		}
	}

	btSequentialImpulseConstraintSolver::solveGroupCacheFriendlySetup(
		bodies, numBodies, manifoldPtr, numManifolds, constraints, numConstraints, infoGlobal, debugDrawer);

//	{
//		BT_PROFILE("gather constraint data");

//		const int numFrictionPerContact = m_tmpSolverContactConstraintPool.size() == m_tmpSolverContactFrictionConstraintPool.size()? 1 : 2;

//	//	int numBodies = m_tmpSolverBodyPool.size();
//		m_allConstraintPtrArray.resize(0);
//		// TODO(JS): Generalize to other constraints
//		m_limitDependencies.resize(m_tmpSolverNonContactConstraintPool.size() + m_tmpSolverContactConstraintPool.size() + m_tmpSolverContactFrictionConstraintPool.size());
////		m_limitDependencies.resize(m_tmpSolverContactConstraintPool.size());
//		btAssert(m_limitDependencies.size() == m_tmpSolverNonContactConstraintPool.size()+m_tmpSolverContactConstraintPool.size()+m_tmpSolverContactFrictionConstraintPool.size());
////		btAssert(m_limitDependencies.size() == m_tmpSolverContactConstraintPool.size());
//	//	printf("m_limitDependencies.size() = %d\n",m_limitDependencies.size());

//		if (m_tmpSolverNonContactConstraintPool.size() != 0)
//			int a = 10;

//		if (m_tmpSolverContactConstraintPool.size() != 0)
//			int a = 10;

//		if (m_tmpSolverContactFrictionConstraintPool.size() != 0)
//			int a = 10;

//		if (m_limitDependencies.size() != 0)
//			int a = 10;

//		int dindex = 0;
//		for (int i = 0; i < m_tmpSolverNonContactConstraintPool.size(); ++i)
//		{
//			m_allConstraintPtrArray.push_back(&m_tmpSolverNonContactConstraintPool[i]);
//			m_limitDependencies[dindex++] = -1;
//		}

//		///The btSequentialImpulseConstraintSolver moves all friction constraints at the very end, we can also interleave them instead

//		int firstContactConstraintOffset=dindex;

//		if (interleaveContactAndFriction)
//		{
//			for (int i = 0; i < m_tmpSolverContactConstraintPool.size(); ++i)
//			{
//				m_allConstraintPtrArray.push_back(&m_tmpSolverContactConstraintPool[i]);
//				m_limitDependencies[dindex++] = -1;
//				m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i*numFrictionPerContact]);
//				const int findex = (m_tmpSolverContactFrictionConstraintPool[i*numFrictionPerContact].m_frictionIndex*(1+numFrictionPerContact));
//				m_limitDependencies[dindex++] = findex + firstContactConstraintOffset;
//				if (numFrictionPerContact == 2)
//				{
//					m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i*numFrictionPerContact+1]);
//					m_limitDependencies[dindex++] = findex + firstContactConstraintOffset;
//				}
//			}
//		}
//		else
//		{
//			for (int i = 0; i < m_tmpSolverContactConstraintPool.size(); ++i)
//			{
//				m_allConstraintPtrArray.push_back(&m_tmpSolverContactConstraintPool[i]);
//				m_limitDependencies[dindex++] = -1;
//			}

//			for (int i = 0; i < m_tmpSolverContactFrictionConstraintPool.size(); ++i)
//			{
//				m_allConstraintPtrArray.push_back(&m_tmpSolverContactFrictionConstraintPool[i]);
//				m_limitDependencies[dindex++] = m_tmpSolverContactFrictionConstraintPool[i].m_frictionIndex + firstContactConstraintOffset;
//			}
//		}

//		if (!m_allConstraintPtrArray.size())
//		{
//			m_A.resize(0, 0);
//			m_b.resize(0);
//			m_x.resize(0);
//			m_lo.resize(0);
//			m_hi.resize(0);

//			return btScalar(0);
//		}
//	}

	{
		BT_PROFILE("createMLCPFast");
		createMLCPFast(infoGlobal);
	}

	return 0;
	// TODO(JS): Not sure what should be returned
}

btScalar btMultiBodyMLCPConstraintSolver::solveGroupCacheFriendlyIterations(btCollisionObject **bodies, int numBodies, btPersistentManifold **manifoldPtr, int numManifolds, btTypedConstraint **constraints, int numConstraints, const btContactSolverInfo &infoGlobal, btIDebugDraw *debugDrawer)
{
	bool success = true;
	{
		BT_PROFILE("solveMLCP");
		printf("m_A(%d,%d)\n", m_A.rows(),m_A.cols());
		print(m_A);
		success = solveMLCP(infoGlobal);
	}

	//check if solution is valid, and otherwise fallback to btSequentialImpulseConstraintSolver::solveGroupCacheFriendlyIterations
	if (success)
	{
		BT_PROFILE("process MLCP results");
		for (int i=0;i<m_allConstraintPtrArray.size();i++)
		{
			{
				btSolverConstraint& c = *m_allConstraintPtrArray[i];
				int sbA = c.m_solverBodyIdA;
				int sbB = c.m_solverBodyIdB;
				//btRigidBody* orgBodyA = m_tmpSolverBodyPool[sbA].m_originalBody;
			//	btRigidBody* orgBodyB = m_tmpSolverBodyPool[sbB].m_originalBody;

				btSolverBody& solverBodyA = m_tmpSolverBodyPool[sbA];
				btSolverBody& solverBodyB = m_tmpSolverBodyPool[sbB];

				{
					btScalar deltaImpulse = m_x[i]-c.m_appliedImpulse;
					c.m_appliedImpulse = m_x[i];
					solverBodyA.internalApplyImpulse(c.m_contactNormal1*solverBodyA.internalGetInvMass(),c.m_angularComponentA,deltaImpulse);
					solverBodyB.internalApplyImpulse(c.m_contactNormal2*solverBodyB.internalGetInvMass(),c.m_angularComponentB,deltaImpulse);
				}

				if (infoGlobal.m_splitImpulse)
				{
					btScalar deltaImpulse = m_xSplit[i] - c.m_appliedPushImpulse;
					solverBodyA.internalApplyPushImpulse(c.m_contactNormal1*solverBodyA.internalGetInvMass(),c.m_angularComponentA,deltaImpulse);
					solverBodyB.internalApplyPushImpulse(c.m_contactNormal2*solverBodyB.internalGetInvMass(),c.m_angularComponentB,deltaImpulse);
					c.m_appliedPushImpulse = m_xSplit[i];
				}

			}
		}
	}
	else
	{
		printf("m_fallback = %d\n",m_fallback);
		m_fallback++;
		btSequentialImpulseConstraintSolver::solveGroupCacheFriendlyIterations(bodies ,numBodies,manifoldPtr, numManifolds,constraints,numConstraints,infoGlobal,debugDrawer);
	}

	return 0.f;
}

void btMultiBodyMLCPConstraintSolver::convertContacts(btPersistentManifold** manifoldPtr, int numManifolds,
													  const btContactSolverInfo& infoGlobal)
{
	for (int i = 0; i < numManifolds; ++i)
	{
		btPersistentManifold* manifold = manifoldPtr[i];
		const btMultiBodyLinkCollider* fcA = btMultiBodyLinkCollider::upcast(manifold->getBody0());
		const btMultiBodyLinkCollider* fcB = btMultiBodyLinkCollider::upcast(manifold->getBody1());
		if (!fcA && !fcB)
		{
			// The contact doesn't involve any Featherstone btMultiBody, so deal with the regular
			// btRigidBody/btCollisionObject case
			convertContact(manifold, infoGlobal);
		}
		else
		{
			convertMultiBodyContact(manifold, infoGlobal);
		}
	}
}

void btMultiBodyMLCPConstraintSolver::convertMultiBodyContact(btPersistentManifold* manifold,
															  const btContactSolverInfo& infoGlobal)
{
	const btMultiBodyLinkCollider* fcA = btMultiBodyLinkCollider::upcast(manifold->getBody0());
	const btMultiBodyLinkCollider* fcB = btMultiBodyLinkCollider::upcast(manifold->getBody1());

	btMultiBody* mbA = fcA ? fcA->m_multiBody :  0;
	btMultiBody* mbB = fcB ? fcB->m_multiBody :  0;

	btCollisionObject* colObj0 =  0;
	btCollisionObject* colObj1 =  0;

	colObj0 = const_cast<btCollisionObject*>(manifold->getBody0());
	colObj1 = const_cast<btCollisionObject*>(manifold->getBody1());

	int solverBodyIdA = mbA ? -1 : getOrInitSolverBody(*colObj0, infoGlobal.m_timeStep);
	int solverBodyIdB = mbB ? -1 : getOrInitSolverBody(*colObj1, infoGlobal.m_timeStep);

	for (int j = 0; j < manifold->getNumContacts(); ++j)
	{
		btManifoldPoint& cp = manifold->getContactPoint(j);

		if (cp.getDistance() <= manifold->getContactProcessingThreshold())
		{
			btScalar relaxation;

			btMultiBodySolverConstraint& solverConstraint = m_multiBodyNormalContactConstraints.expandNonInitializing();
			solverConstraint.m_orgConstraint = 0;
			solverConstraint.m_orgDofIndex = -1;
			solverConstraint.m_solverBodyIdA = solverBodyIdA;
			solverConstraint.m_solverBodyIdB = solverBodyIdB;
			solverConstraint.m_multiBodyA = mbA;
			if (mbA) solverConstraint.m_linkA = fcA->m_link;
			solverConstraint.m_multiBodyB = mbB;
			if (mbB) solverConstraint.m_linkB = fcB->m_link;
			solverConstraint.m_originalContactPoint = &cp;

			const bool isFriction = false;
			setupMultiBodyContactConstraint(solverConstraint, cp.m_normalWorldOnB, cp, infoGlobal, relaxation,
											isFriction);
		}
	}
}

void btMultiBodyMLCPConstraintSolver::setupMultiBodyContactConstraint(
	btMultiBodySolverConstraint& solverConstraint, const btVector3& contactNormal, btManifoldPoint& cp,
	const btContactSolverInfo& infoGlobal, btScalar& relaxation, bool isFriction, btScalar /*desiredVelocity*/,
	btScalar /*cfmSlip*/)
{
	BT_PROFILE("btMultiBodyMLCPConstraintSolver::setupMultiBodyContactConstraint");

	btMultiBody* multiBodyA = solverConstraint.m_multiBodyA;
	btMultiBody* multiBodyB = solverConstraint.m_multiBodyB;

	const btVector3& contactPositionOnA = cp.getPositionWorldOnA();
	const btVector3& contactPositionOnB = cp.getPositionWorldOnB();

	btSolverBody* bodyA = multiBodyA ?  0 : &m_tmpSolverBodyPool[solverConstraint.m_solverBodyIdA];
	btSolverBody* bodyB = multiBodyB ?  0 : &m_tmpSolverBodyPool[solverConstraint.m_solverBodyIdB];

	btRigidBody* rbA = multiBodyA ?  0 : bodyA->m_originalBody;
	btRigidBody* rbB = multiBodyB ?  0 : bodyB->m_originalBody;

	// Relative position from the origin of bodyA to the contact point on bodyA
	btVector3 relPosA;
	if (bodyA) relPosA = contactPositionOnA - bodyA->getWorldTransform().getOrigin();

	// Relative position from the origin of bodyB to the contact point on bodyB
	btVector3 relPosB;
	if (bodyB) relPosB = contactPositionOnB - bodyB->getWorldTransform().getOrigin();

	relaxation = infoGlobal.m_sor;

	const btScalar invTimeStep = btScalar(1) / infoGlobal.m_timeStep;

	btScalar cfm;
	btScalar erp;
	if (isFriction)
	{
		cfm = infoGlobal.m_frictionCFM;
		erp = infoGlobal.m_frictionERP;
	}
	else
	{
		cfm = infoGlobal.m_globalCfm;
		erp = infoGlobal.m_erp2;

		if ((cp.m_contactPointFlags & BT_CONTACT_FLAG_HAS_CONTACT_CFM) ||
			(cp.m_contactPointFlags & BT_CONTACT_FLAG_HAS_CONTACT_ERP))
		{
			if (cp.m_contactPointFlags & BT_CONTACT_FLAG_HAS_CONTACT_CFM) cfm = cp.m_contactCFM;
			if (cp.m_contactPointFlags & BT_CONTACT_FLAG_HAS_CONTACT_ERP) erp = cp.m_contactERP;
		}
		else
		{
			if (cp.m_contactPointFlags & BT_CONTACT_FLAG_CONTACT_STIFFNESS_DAMPING)
			{
				btScalar denom =
					(infoGlobal.m_timeStep * cp.m_combinedContactStiffness1 + cp.m_combinedContactDamping1);
				if (denom < SIMD_EPSILON) denom = SIMD_EPSILON;
				cfm = btScalar(1) / denom;
				erp = (infoGlobal.m_timeStep * cp.m_combinedContactStiffness1) / denom;
			}
		}
	}

	cfm *= invTimeStep;

	if (multiBodyA)
	{
		if (solverConstraint.m_linkA < 0)
		{
			relPosA = contactPositionOnA - multiBodyA->getBasePos();
		}
		else
		{
			relPosA =
				contactPositionOnA - multiBodyA->getLink(solverConstraint.m_linkA).m_cachedWorldTransform.getOrigin();
			// TODO(JS): Does this can be different from bodyA->getWorldTransform().getOrigin()?
		}
		const int ndofA = multiBodyA->getNumDofs() + 6;

		// Append joint velocity change to the array, if not already, and set the index to it for MultiBody A
		solverConstraint.m_deltaVelAindex = multiBodyA->getCompanionId();
		if (solverConstraint.m_deltaVelAindex < 0)
		{
			solverConstraint.m_deltaVelAindex = m_data.m_deltaVelocities.size();
			multiBodyA->setCompanionId(solverConstraint.m_deltaVelAindex);
			m_data.m_deltaVelocities.resize(m_data.m_deltaVelocities.size() + ndofA);
		}
		else
		{
			btAssert(m_data.m_deltaVelocities.size() >= solverConstraint.m_deltaVelAindex + ndofA);
		}

		// Append contact Jacobians to the array, if not already, and set the index to it for MultiBody A
		solverConstraint.m_jacAindex = m_data.m_jacobians.size();
		m_data.m_jacobians.resize(m_data.m_jacobians.size() + ndofA);
		m_data.m_deltaVelocitiesUnitImpulse.resize(m_data.m_deltaVelocitiesUnitImpulse.size() + ndofA);
		btAssert(m_data.m_jacobians.size() == m_data.m_deltaVelocitiesUnitImpulse.size());

		btScalar* jac1 = &m_data.m_jacobians[solverConstraint.m_jacAindex];
		multiBodyA->fillContactJacobianMultiDof(solverConstraint.m_linkA, cp.getPositionWorldOnA(), contactNormal, jac1,
												m_data.scratch_r, m_data.scratch_v, m_data.scratch_m);
		btScalar* delta = &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacAindex];
		multiBodyA->calcAccelerationDeltasMultiDof(&m_data.m_jacobians[solverConstraint.m_jacAindex], delta,
												   m_data.scratch_r, m_data.scratch_v);
		// NOTE(JS): The velocity changes in other constraint spaces are necessary to be stored for MultibodyMLCP solver
		//           So this should be moved to btMultiBodyMLCPSolver::createMultiBodyMCLPFast()

		const btVector3 torqueAxis0 = relPosA.cross(contactNormal);
		solverConstraint.m_relpos1CrossNormal = torqueAxis0;
		solverConstraint.m_contactNormal1 = contactNormal;
	}
	else
	{
		const btVector3 torqueAxis0 = relPosA.cross(contactNormal);
		solverConstraint.m_relpos1CrossNormal = torqueAxis0;
		solverConstraint.m_contactNormal1 = contactNormal;
		solverConstraint.m_angularComponentA =
			rbA ? rbA->getInvInertiaTensorWorld() * torqueAxis0 * rbA->getAngularFactor() : btVector3(0, 0, 0);
	}

	if (multiBodyB)
	{
		if (solverConstraint.m_linkB < 0)
		{
			relPosB = contactPositionOnB - multiBodyB->getBasePos();
		}
		else
		{
			relPosB =
				contactPositionOnB - multiBodyB->getLink(solverConstraint.m_linkB).m_cachedWorldTransform.getOrigin();
		}

		const int ndofB = multiBodyB->getNumDofs() + 6;

		solverConstraint.m_deltaVelBindex = multiBodyB->getCompanionId();
		if (solverConstraint.m_deltaVelBindex < 0)
		{
			solverConstraint.m_deltaVelBindex = m_data.m_deltaVelocities.size();
			multiBodyB->setCompanionId(solverConstraint.m_deltaVelBindex);
			m_data.m_deltaVelocities.resize(m_data.m_deltaVelocities.size() + ndofB);
		}

		solverConstraint.m_jacBindex = m_data.m_jacobians.size();

		m_data.m_jacobians.resize(m_data.m_jacobians.size() + ndofB);
		m_data.m_deltaVelocitiesUnitImpulse.resize(m_data.m_deltaVelocitiesUnitImpulse.size() + ndofB);
		btAssert(m_data.m_jacobians.size() == m_data.m_deltaVelocitiesUnitImpulse.size());

		multiBodyB->fillContactJacobianMultiDof(solverConstraint.m_linkB, cp.getPositionWorldOnB(), -contactNormal,
												&m_data.m_jacobians[solverConstraint.m_jacBindex], m_data.scratch_r,
												m_data.scratch_v, m_data.scratch_m);
		multiBodyB->calcAccelerationDeltasMultiDof(&m_data.m_jacobians[solverConstraint.m_jacBindex],
												   &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacBindex],
												   m_data.scratch_r, m_data.scratch_v);

		const btVector3 torqueAxis1 = relPosB.cross(contactNormal);
		solverConstraint.m_relpos2CrossNormal = -torqueAxis1;
		solverConstraint.m_contactNormal2 = -contactNormal;
	}
	else
	{
		const btVector3 torqueAxis1 = relPosB.cross(contactNormal);
		solverConstraint.m_relpos2CrossNormal = -torqueAxis1;
		solverConstraint.m_contactNormal2 = -contactNormal;

		solverConstraint.m_angularComponentB =
			rbB ? rbB->getInvInertiaTensorWorld() * -torqueAxis1 * rbB->getAngularFactor() : btVector3(0, 0, 0);
	}

	{
		btVector3 vec;
		btScalar denom0 = btScalar(0);
		btScalar denom1 = btScalar(0);
		btScalar* jacB =  0;
		btScalar* jacA =  0;
		btScalar* lambdaA =  0;
		btScalar* lambdaB =  0;
		int ndofA = 0;
		if (multiBodyA)
		{
			ndofA = multiBodyA->getNumDofs() + 6;
			jacA = &m_data.m_jacobians[solverConstraint.m_jacAindex];
			lambdaA = &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacAindex];
			for (int i = 0; i < ndofA; ++i)
			{
				const btScalar j = jacA[i];
				const btScalar l = lambdaA[i];
				denom0 += j * l;
			}
		}
		else
		{
			if (rbA)
			{
				vec = (solverConstraint.m_angularComponentA).cross(relPosA);
				denom0 = rbA->getInvMass() + contactNormal.dot(vec);
			}
		}
		if (multiBodyB)
		{
			const int ndofB = multiBodyB->getNumDofs() + 6;
			jacB = &m_data.m_jacobians[solverConstraint.m_jacBindex];
			lambdaB = &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacBindex];
			for (int i = 0; i < ndofB; ++i)
			{
				const btScalar j = jacB[i];
				const btScalar l = lambdaB[i];
				denom1 += j * l;
			}
		}
		else
		{
			if (rbB)
			{
				vec = (-solverConstraint.m_angularComponentB).cross(relPosB);
				denom1 = rbB->getInvMass() + contactNormal.dot(vec);
			}
		}

		const btScalar d = denom0 + denom1 + cfm;
		if (d > SIMD_EPSILON)
		{
			solverConstraint.m_jacDiagABInv = relaxation / d;
		}
		else
		{
			//disable the constraint row to handle singularity/redundant constraint
			solverConstraint.m_jacDiagABInv = btScalar(0);
		}
	}

	//compute rhs and remaining solverConstraint fields

	btScalar restitution = btScalar(0);
	btScalar distance = btScalar(0);
	if (!isFriction)
	{
		distance = cp.getDistance() + infoGlobal.m_linearSlop;
	}
	else
	{
		if (cp.m_contactPointFlags & BT_CONTACT_FLAG_FRICTION_ANCHOR)
		{
			distance = (cp.getPositionWorldOnA() - cp.getPositionWorldOnB()).dot(contactNormal);
		}
	}

	btScalar rel_vel = btScalar(0);
	int ndofA = 0;
	int ndofB = 0;
	{
		if (multiBodyA)
		{
			ndofA = multiBodyA->getNumDofs() + 6;
			btScalar* jacA = &m_data.m_jacobians[solverConstraint.m_jacAindex];
			for (int i = 0; i < ndofA; ++i) rel_vel += multiBodyA->getVelocityVector()[i] * jacA[i];
		}
		else
		{
			if (rbA)
			{
				rel_vel +=
					(rbA->getVelocityInLocalPoint(relPosA) +
					 (rbA->getTotalTorque() * rbA->getInvInertiaTensorWorld() * infoGlobal.m_timeStep).cross(relPosA) +
					 rbA->getTotalForce() * rbA->getInvMass() * infoGlobal.m_timeStep)
						.dot(solverConstraint.m_contactNormal1);
			}
		}
		if (multiBodyB)
		{
			ndofB = multiBodyB->getNumDofs() + 6;
			btScalar* jacB = &m_data.m_jacobians[solverConstraint.m_jacBindex];
			for (int i = 0; i < ndofB; ++i) rel_vel += multiBodyB->getVelocityVector()[i] * jacB[i];
		}
		else
		{
			if (rbB)
			{
				rel_vel +=
					(rbB->getVelocityInLocalPoint(relPosB) +
					 (rbB->getTotalTorque() * rbB->getInvInertiaTensorWorld() * infoGlobal.m_timeStep).cross(relPosB) +
					 rbB->getTotalForce() * rbB->getInvMass() * infoGlobal.m_timeStep)
						.dot(solverConstraint.m_contactNormal2);
			}
		}

		solverConstraint.m_friction = cp.m_combinedFriction;

		if (!isFriction)
		{
			restitution =
				restitutionCurve(rel_vel, cp.m_combinedRestitution, infoGlobal.m_restitutionVelocityThreshold);
			restitution = std::max(btScalar(0), restitution);
		}
	}

	///warm starting (or zero if disabled)
	//disable warmstarting for btMultiBody, it has issues gaining energy (==explosion)
	if (0)  //infoGlobal.m_solverMode & SOLVER_USE_WARMSTARTING)
	{
		solverConstraint.m_appliedImpulse = isFriction ? 0 : cp.m_appliedImpulse * infoGlobal.m_warmstartingFactor;

		if (solverConstraint.m_appliedImpulse)
		{
			if (multiBodyA)
			{
				const btScalar impulse = solverConstraint.m_appliedImpulse;
				btScalar* deltaV = &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacAindex];
				multiBodyA->applyDeltaVeeMultiDof(deltaV, impulse);

				applyDeltaVee(deltaV, impulse, solverConstraint.m_deltaVelAindex, ndofA);
			}
			else
			{
				if (rbA)
					bodyA->internalApplyImpulse(
						solverConstraint.m_contactNormal1 * bodyA->internalGetInvMass() * rbA->getLinearFactor(),
						solverConstraint.m_angularComponentA, solverConstraint.m_appliedImpulse);
			}
			if (multiBodyB)
			{
				const btScalar impulse = solverConstraint.m_appliedImpulse;
				btScalar* deltaV = &m_data.m_deltaVelocitiesUnitImpulse[solverConstraint.m_jacBindex];
				multiBodyB->applyDeltaVeeMultiDof(deltaV, impulse);
				applyDeltaVee(deltaV, impulse, solverConstraint.m_deltaVelBindex, ndofB);
			}
			else
			{
				if (rbB)
					bodyB->internalApplyImpulse(
						-solverConstraint.m_contactNormal2 * bodyB->internalGetInvMass() * rbB->getLinearFactor(),
						-solverConstraint.m_angularComponentB, -(btScalar)solverConstraint.m_appliedImpulse);
			}
		}
	}
	else
	{
		solverConstraint.m_appliedImpulse = 0.f;
	}

	solverConstraint.m_appliedPushImpulse = 0.f;

	{
		btScalar positionalError = btScalar(0);
		btScalar velocityError =
			restitution -
			rel_vel;  // * damping;	//note for friction restitution is always set to 0 (check above) so it is acutally velocityError = -rel_vel for friction
		if (isFriction)
		{
			positionalError = -distance * erp / infoGlobal.m_timeStep;
		}
		else
		{
			if (distance > 0)
			{
				positionalError = 0;
				velocityError -= distance / infoGlobal.m_timeStep;
			}
			else
			{
				positionalError = -distance * erp / infoGlobal.m_timeStep;
			}
		}

		const btScalar penetrationImpulse = positionalError * solverConstraint.m_jacDiagABInv;
		const btScalar velocityImpulse = velocityError * solverConstraint.m_jacDiagABInv;

		if (!isFriction)
		{
			//	if (!infoGlobal.m_splitImpulse || (penetration > infoGlobal.m_splitImpulsePenetrationThreshold))
			{
				//combine position and velocity into rhs
				solverConstraint.m_rhs = penetrationImpulse + velocityImpulse;
				solverConstraint.m_rhsPenetration = 0.f;
			}
			/*else
			{
				//split position and velocity into rhs and m_rhsPenetration
				solverConstraint.m_rhs = velocityImpulse;
				solverConstraint.m_rhsPenetration = penetrationImpulse;
			}
			*/
			solverConstraint.m_lowerLimit = 0;
			solverConstraint.m_upperLimit = 1e10f;
		}
		else
		{
			solverConstraint.m_rhs = penetrationImpulse + velocityImpulse;
			solverConstraint.m_rhsPenetration = 0.f;
			solverConstraint.m_lowerLimit = -solverConstraint.m_friction;
			solverConstraint.m_upperLimit = solverConstraint.m_friction;
		}

		solverConstraint.m_cfm = cfm * solverConstraint.m_jacDiagABInv;
	}
}

void btMultiBodyMLCPConstraintSolver::applyDeltaVee(btScalar* deltaV, btScalar impulse, int velocityIndex, int ndof)
{
	// TODO(JS):
}

btMultiBodyMLCPConstraintSolver::btMultiBodyMLCPConstraintSolver(btMLCPSolverInterface *solver)
	: m_solver(solver)
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
