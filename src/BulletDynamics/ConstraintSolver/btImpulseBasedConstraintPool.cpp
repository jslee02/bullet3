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

//#define COMPUTE_IMPULSE_DENOM 1
//#define BT_ADDITIONAL_DEBUG

//It is not necessary (redundant) to refresh contact manifolds, this refresh has been moved to the collision algorithms.

#include "btImpulseBasedConstraintPool.h"

#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"

btScalar btImpulseBasedConstraintPool::setup(btCollisionObject** bodies, int numBodies, btPersistentManifold** manifoldPtr, int numManifolds, btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal, btIDebugDraw* debugDrawer)
{
	//convert all bodies
	convertBodies(bodies, numBodies, infoGlobal);

	convertJoints(constraints, numConstraints, infoGlobal);

	convertContacts(manifoldPtr, numManifolds, infoGlobal);

	//	btContactSolverInfo info = infoGlobal;

	int numNonContactPool = m_tmpSolverNonContactConstraintPool.size();
	int numConstraintPool = m_tmpSolverContactConstraintPool.size();
	int numFrictionPool = m_tmpSolverContactFrictionConstraintPool.size();

	///@todo: use stack allocator for such temporarily memory, same for solver bodies/constraints
	m_orderNonContactConstraintPool.resizeNoInitialize(numNonContactPool);
	if ((infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS))
		m_orderTmpConstraintPool.resizeNoInitialize(numConstraintPool * 2);
	else
		m_orderTmpConstraintPool.resizeNoInitialize(numConstraintPool);

	m_orderFrictionConstraintPool.resizeNoInitialize(numFrictionPool);
	{
		int i;
		for (i = 0; i < numNonContactPool; i++)
		{
			m_orderNonContactConstraintPool[i] = i;
		}
		for (i = 0; i < numConstraintPool; i++)
		{
			m_orderTmpConstraintPool[i] = i;
		}
		for (i = 0; i < numFrictionPool; i++)
		{
			m_orderFrictionConstraintPool[i] = i;
		}
	}

	return 0.f;
}

void btImpulseBasedConstraintPool::convertContacts(btPersistentManifold** manifoldPtr, int numManifolds, const btContactSolverInfo& infoGlobal)
{
	int i;
	btPersistentManifold* manifold = 0;
	//			btCollisionObject* colObj0=0,*colObj1=0;

	for (i = 0; i < numManifolds; i++)
	{
		manifold = manifoldPtr[i];
		convertContact(manifold, infoGlobal);
	}
}

void btImpulseBasedConstraintPool::convertContact(btPersistentManifold* manifold, const btContactSolverInfo& infoGlobal)
{
	btCollisionObject *colObj0 = 0, *colObj1 = 0;

	colObj0 = (btCollisionObject*)manifold->getBody0();
	colObj1 = (btCollisionObject*)manifold->getBody1();

	int solverBodyIdA = getOrInitSolverBody(*colObj0, infoGlobal.m_timeStep);
	int solverBodyIdB = getOrInitSolverBody(*colObj1, infoGlobal.m_timeStep);

	//	btRigidBody* bodyA = btRigidBody::upcast(colObj0);
	//	btRigidBody* bodyB = btRigidBody::upcast(colObj1);

	btSolverBody* solverBodyA = &m_tmpSolverBodyPool[solverBodyIdA];
	btSolverBody* solverBodyB = &m_tmpSolverBodyPool[solverBodyIdB];

	///avoid collision response between two static objects
	if (!solverBodyA || (solverBodyA->m_invMass.fuzzyZero() && (!solverBodyB || solverBodyB->m_invMass.fuzzyZero())))
		return;

	int rollingFriction = 1;
	for (int j = 0; j < manifold->getNumContacts(); j++)
	{
		btManifoldPoint& cp = manifold->getContactPoint(j);

		if (cp.getDistance() <= manifold->getContactProcessingThreshold())
		{
			btVector3 rel_pos1;
			btVector3 rel_pos2;
			btScalar relaxation;

			int frictionIndex = m_tmpSolverContactConstraintPool.size();
			btSolverConstraint& solverConstraint = m_tmpSolverContactConstraintPool.expandNonInitializing();
			solverConstraint.m_solverBodyIdA = solverBodyIdA;
			solverConstraint.m_solverBodyIdB = solverBodyIdB;

			solverConstraint.m_originalContactPoint = &cp;

			const btVector3& pos1 = cp.getPositionWorldOnA();
			const btVector3& pos2 = cp.getPositionWorldOnB();

			rel_pos1 = pos1 - colObj0->getWorldTransform().getOrigin();
			rel_pos2 = pos2 - colObj1->getWorldTransform().getOrigin();

			btVector3 vel1;
			btVector3 vel2;

			solverBodyA->getVelocityInLocalPointNoDelta(rel_pos1, vel1);
			solverBodyB->getVelocityInLocalPointNoDelta(rel_pos2, vel2);

			btVector3 vel = vel1 - vel2;
			btScalar rel_vel = cp.m_normalWorldOnB.dot(vel);

			setupContactConstraint(solverConstraint, solverBodyIdA, solverBodyIdB, cp, infoGlobal, relaxation, rel_pos1, rel_pos2);

			/////setup the friction constraints

			solverConstraint.m_frictionIndex = m_tmpSolverContactFrictionConstraintPool.size();

			if ((cp.m_combinedRollingFriction > 0.f) && (rollingFriction > 0))
			{
				{
					addTorsionalFrictionConstraint(cp.m_normalWorldOnB, solverBodyIdA, solverBodyIdB, frictionIndex, cp, cp.m_combinedSpinningFriction, rel_pos1, rel_pos2, colObj0, colObj1, relaxation);
					btVector3 axis0, axis1;
					btPlaneSpace1(cp.m_normalWorldOnB, axis0, axis1);
					axis0.normalize();
					axis1.normalize();

					applyAnisotropicFriction(colObj0, axis0, btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION);
					applyAnisotropicFriction(colObj1, axis0, btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION);
					applyAnisotropicFriction(colObj0, axis1, btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION);
					applyAnisotropicFriction(colObj1, axis1, btCollisionObject::CF_ANISOTROPIC_ROLLING_FRICTION);
					if (axis0.length() > 0.001)
						addTorsionalFrictionConstraint(axis0, solverBodyIdA, solverBodyIdB, frictionIndex, cp,
													   cp.m_combinedRollingFriction, rel_pos1, rel_pos2, colObj0, colObj1, relaxation);
					if (axis1.length() > 0.001)
						addTorsionalFrictionConstraint(axis1, solverBodyIdA, solverBodyIdB, frictionIndex, cp,
													   cp.m_combinedRollingFriction, rel_pos1, rel_pos2, colObj0, colObj1, relaxation);
				}
			}

			///Bullet has several options to set the friction directions
			///By default, each contact has only a single friction direction that is recomputed automatically very frame
			///based on the relative linear velocity.
			///If the relative velocity it zero, it will automatically compute a friction direction.

			///You can also enable two friction directions, using the SOLVER_USE_2_FRICTION_DIRECTIONS.
			///In that case, the second friction direction will be orthogonal to both contact normal and first friction direction.
			///
			///If you choose SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION, then the friction will be independent from the relative projected velocity.
			///
			///The user can manually override the friction directions for certain contacts using a contact callback,
			///and set the cp.m_lateralFrictionInitialized to true
			///In that case, you can set the target relative motion in each friction direction (cp.m_contactMotion1 and cp.m_contactMotion2)
			///this will give a conveyor belt effect
			///

			if (!(infoGlobal.m_solverMode & SOLVER_ENABLE_FRICTION_DIRECTION_CACHING) || !(cp.m_contactPointFlags & BT_CONTACT_FLAG_LATERAL_FRICTION_INITIALIZED))
			{
				cp.m_lateralFrictionDir1 = vel - cp.m_normalWorldOnB * rel_vel;
				btScalar lat_rel_vel = cp.m_lateralFrictionDir1.length2();
				if (!(infoGlobal.m_solverMode & SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION) && lat_rel_vel > SIMD_EPSILON)
				{
					cp.m_lateralFrictionDir1 *= 1.f / btSqrt(lat_rel_vel);
					applyAnisotropicFriction(colObj0, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION);
					applyAnisotropicFriction(colObj1, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION);
					addFrictionConstraint(cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal);

					if ((infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS))
					{
						cp.m_lateralFrictionDir2 = cp.m_lateralFrictionDir1.cross(cp.m_normalWorldOnB);
						cp.m_lateralFrictionDir2.normalize();  //??
						applyAnisotropicFriction(colObj0, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION);
						applyAnisotropicFriction(colObj1, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION);
						addFrictionConstraint(cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal);
					}
				}
				else
				{
					btPlaneSpace1(cp.m_normalWorldOnB, cp.m_lateralFrictionDir1, cp.m_lateralFrictionDir2);

					applyAnisotropicFriction(colObj0, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION);
					applyAnisotropicFriction(colObj1, cp.m_lateralFrictionDir1, btCollisionObject::CF_ANISOTROPIC_FRICTION);
					addFrictionConstraint(cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal);

					if ((infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS))
					{
						applyAnisotropicFriction(colObj0, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION);
						applyAnisotropicFriction(colObj1, cp.m_lateralFrictionDir2, btCollisionObject::CF_ANISOTROPIC_FRICTION);
						addFrictionConstraint(cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal);
					}

					if ((infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS) && (infoGlobal.m_solverMode & SOLVER_DISABLE_VELOCITY_DEPENDENT_FRICTION_DIRECTION))
					{
						cp.m_contactPointFlags |= BT_CONTACT_FLAG_LATERAL_FRICTION_INITIALIZED;
					}
				}
			}
			else
			{
				addFrictionConstraint(cp.m_lateralFrictionDir1, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal, cp.m_contactMotion1, cp.m_frictionCFM);

				if ((infoGlobal.m_solverMode & SOLVER_USE_2_FRICTION_DIRECTIONS))
					addFrictionConstraint(cp.m_lateralFrictionDir2, solverBodyIdA, solverBodyIdB, frictionIndex, cp, rel_pos1, rel_pos2, colObj0, colObj1, relaxation, infoGlobal, cp.m_contactMotion2, cp.m_frictionCFM);
			}
			setFrictionConstraintImpulse(solverConstraint, solverBodyIdA, solverBodyIdB, cp, infoGlobal);
		}
	}
}

btSolverConstraint &btImpulseBasedConstraintPool::addFrictionConstraint(const btVector3 &normalAxis, int solverBodyIdA, int solverBodyIdB, int frictionIndex, btManifoldPoint &cp, const btVector3 &rel_pos1, const btVector3 &rel_pos2, btCollisionObject *colObj0, btCollisionObject *colObj1, btScalar relaxation, const btContactSolverInfo &infoGlobal, btScalar desiredVelocity, btScalar cfmSlip)
{

}

btSolverConstraint &btImpulseBasedConstraintPool::addTorsionalFrictionConstraint(const btVector3 &normalAxis, int solverBodyIdA, int solverBodyIdB, int frictionIndex, btManifoldPoint &cp, btScalar torsionalFriction, const btVector3 &rel_pos1, const btVector3 &rel_pos2, btCollisionObject *colObj0, btCollisionObject *colObj1, btScalar relaxation, btScalar desiredVelocity, btScalar cfmSlip)
{

}

void btImpulseBasedConstraintPool::applyAnisotropicFriction(btCollisionObject* colObj, btVector3& frictionDirection, int frictionMode)
{
	if (colObj && colObj->hasAnisotropicFriction(frictionMode))
	{
		// transform to local coordinates
		btVector3 loc_lateral = frictionDirection * colObj->getWorldTransform().getBasis();
		const btVector3& friction_scaling = colObj->getAnisotropicFriction();
		//apply anisotropic friction
		loc_lateral *= friction_scaling;
		// ... and transform it back to global coordinates
		frictionDirection = colObj->getWorldTransform().getBasis() * loc_lateral;
	}
}

void btImpulseBasedConstraintPool::convertJoints(btTypedConstraint** constraints, int numConstraints, const btContactSolverInfo& infoGlobal)
{
}

void btImpulseBasedConstraintPool::convertJoint(btSolverConstraint* currentConstraintRow, btTypedConstraint* constraint, const btTypedConstraint::btConstraintInfo1& info1, int solverBodyIdA, int solverBodyIdB, const btContactSolverInfo& infoGlobal)
{
}

void btImpulseBasedConstraintPool::convertBodies(btCollisionObject** bodies, int numBodies, const btContactSolverInfo& infoGlobal)
{
}

int btImpulseBasedConstraintPool::getOrInitSolverBody(btCollisionObject& body, btScalar timeStep)
{
#if BT_THREADSAFE
	int solverBodyId = -1;
	bool isRigidBodyType = btRigidBody::upcast(&body) != NULL;
	if (isRigidBodyType && !body.isStaticOrKinematicObject())
	{
		// dynamic body
		// Dynamic bodies can only be in one island, so it's safe to write to the companionId
		solverBodyId = body.getCompanionId();
		if (solverBodyId < 0)
		{
			solverBodyId = m_tmpSolverBodyPool.size();
			btSolverBody& solverBody = m_tmpSolverBodyPool.expand();
			initSolverBody(&solverBody, &body, timeStep);
			body.setCompanionId(solverBodyId);
		}
	}
	else if (isRigidBodyType && body.isKinematicObject())
	{
		//
		// NOTE: must test for kinematic before static because some kinematic objects also
		//   identify as "static"
		//
		// Kinematic bodies can be in multiple islands at once, so it is a
		// race condition to write to them, so we use an alternate method
		// to record the solverBodyId
		int uniqueId = body.getWorldArrayIndex();
		const int INVALID_SOLVER_BODY_ID = -1;
		if (uniqueId >= m_kinematicBodyUniqueIdToSolverBodyTable.size())
		{
			m_kinematicBodyUniqueIdToSolverBodyTable.resize(uniqueId + 1, INVALID_SOLVER_BODY_ID);
		}
		solverBodyId = m_kinematicBodyUniqueIdToSolverBodyTable[uniqueId];
		// if no table entry yet,
		if (solverBodyId == INVALID_SOLVER_BODY_ID)
		{
			// create a table entry for this body
			solverBodyId = m_tmpSolverBodyPool.size();
			btSolverBody& solverBody = m_tmpSolverBodyPool.expand();
			initSolverBody(&solverBody, &body, timeStep);
			m_kinematicBodyUniqueIdToSolverBodyTable[uniqueId] = solverBodyId;
		}
	}
	else
	{
		bool isMultiBodyType = (body.getInternalType() & btCollisionObject::CO_FEATHERSTONE_LINK);
		// Incorrectly set collision object flags can degrade performance in various ways.
		if (!isMultiBodyType)
		{
			btAssert(body.isStaticOrKinematicObject());
		}
		//it could be a multibody link collider
		// all fixed bodies (inf mass) get mapped to a single solver id
		if (m_fixedBodyId < 0)
		{
			m_fixedBodyId = m_tmpSolverBodyPool.size();
			btSolverBody& fixedBody = m_tmpSolverBodyPool.expand();
			initSolverBody(&fixedBody, 0, timeStep);
		}
		solverBodyId = m_fixedBodyId;
	}
	btAssert(solverBodyId >= 0 && solverBodyId < m_tmpSolverBodyPool.size());
	return solverBodyId;
#else   // BT_THREADSAFE

	int solverBodyIdA = -1;

	if (body.getCompanionId() >= 0)
	{
		//body has already been converted
		solverBodyIdA = body.getCompanionId();
		btAssert(solverBodyIdA < m_tmpSolverBodyPool.size());
	}
	else
	{
		btRigidBody* rb = btRigidBody::upcast(&body);
		//convert both active and kinematic objects (for their velocity)
		if (rb && (rb->getInvMass() || rb->isKinematicObject()))
		{
			solverBodyIdA = m_tmpSolverBodyPool.size();
			btSolverBody& solverBody = m_tmpSolverBodyPool.expand();
			initSolverBody(&solverBody, &body, timeStep);
			body.setCompanionId(solverBodyIdA);
		}
		else
		{
			if (m_fixedBodyId < 0)
			{
				m_fixedBodyId = m_tmpSolverBodyPool.size();
				btSolverBody& fixedBody = m_tmpSolverBodyPool.expand();
				initSolverBody(&fixedBody, 0, timeStep);
			}
			return m_fixedBodyId;
			//			return 0;//assume first one is a fixed solver body
		}
	}

	return solverBodyIdA;
#endif  // BT_THREADSAFE
}

void btImpulseBasedConstraintPool::initSolverBody(btSolverBody* solverBody, btCollisionObject* collisionObject, btScalar timeStep)
{
	btRigidBody* rb = collisionObject ? btRigidBody::upcast(collisionObject) : 0;

	solverBody->internalGetDeltaLinearVelocity().setValue(0.f, 0.f, 0.f);
	solverBody->internalGetDeltaAngularVelocity().setValue(0.f, 0.f, 0.f);
	solverBody->internalGetPushVelocity().setValue(0.f, 0.f, 0.f);
	solverBody->internalGetTurnVelocity().setValue(0.f, 0.f, 0.f);

	if (rb)
	{
		solverBody->m_worldTransform = rb->getWorldTransform();
		solverBody->internalSetInvMass(btVector3(rb->getInvMass(), rb->getInvMass(), rb->getInvMass()) * rb->getLinearFactor());
		solverBody->m_originalBody = rb;
		solverBody->m_angularFactor = rb->getAngularFactor();
		solverBody->m_linearFactor = rb->getLinearFactor();
		solverBody->m_linearVelocity = rb->getLinearVelocity();
		solverBody->m_angularVelocity = rb->getAngularVelocity();
		solverBody->m_externalForceImpulse = rb->getTotalForce() * rb->getInvMass() * timeStep;
		solverBody->m_externalTorqueImpulse = rb->getTotalTorque() * rb->getInvInertiaTensorWorld() * timeStep;
	}
	else
	{
		solverBody->m_worldTransform.setIdentity();
		solverBody->internalSetInvMass(btVector3(0, 0, 0));
		solverBody->m_originalBody = 0;
		solverBody->m_angularFactor.setValue(1, 1, 1);
		solverBody->m_linearFactor.setValue(1, 1, 1);
		solverBody->m_linearVelocity.setValue(0, 0, 0);
		solverBody->m_angularVelocity.setValue(0, 0, 0);
		solverBody->m_externalForceImpulse.setValue(0, 0, 0);
		solverBody->m_externalTorqueImpulse.setValue(0, 0, 0);
	}
}
