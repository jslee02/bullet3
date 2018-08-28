/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2013 Erwin Coumans  http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef BT_MULTIBODY_CONSTRAINT_SOLVER_H
#define BT_MULTIBODY_CONSTRAINT_SOLVER_H

#include "BulletDynamics/ConstraintSolver/btSequentialImpulseConstraintSolver.h"
#include "btMultiBodySolverConstraint.h"

#define DIRECTLY_UPDATE_VELOCITY_DURING_SOLVER_ITERATIONS

class btMultiBody;

#include "btMultiBodyConstraint.h"



ATTRIBUTE_ALIGNED16(class) btMultiBodyConstraintSolver : public btSequentialImpulseConstraintSolver
{

protected:

	btMultiBodyConstraintArray			m_multiBodyNonContactConstraints;

	btMultiBodyConstraintArray			m_multiBodyNormalContactConstraints;
	btMultiBodyConstraintArray			m_multiBodyFrictionContactConstraints;
	btMultiBodyConstraintArray			m_multiBodyTorsionalFrictionContactConstraints;

	btMultiBodyJacobianData				m_data;
	
	//temp storage for multi body constraints for a specific island/group called by 'solveGroup'
	btMultiBodyConstraint**					m_tmpMultiBodyConstraints;
	int										m_tmpNumMultiBodyConstraints;

	btScalar resolveSingleConstraintRowGeneric(const btMultiBodySolverConstraint& c);
	
	//solve 2 friction directions and clamp against the implicit friction cone
	btScalar resolveConeFrictionConstraintRows(const btMultiBodySolverConstraint& cA1, const btMultiBodySolverConstraint& cB);
	

	void convertContacts(btPersistentManifold** manifoldPtr,int numManifolds, const btContactSolverInfo& infoGlobal);
    
    btMultiBodySolverConstraint&	addMultiBodyFrictionConstraint(const btVector3& normalAxis,btPersistentManifold* manifold,int frictionIndex,btManifoldPoint& cp,btCollisionObject* colObj0,btCollisionObject* colObj1, btScalar relaxation, const btContactSolverInfo& infoGlobal, btScalar desiredVelocity=0, btScalar cfmSlip=0);

    btMultiBodySolverConstraint&	addMultiBodyTorsionalFrictionConstraint(const btVector3& normalAxis,btPersistentManifold* manifold,int frictionIndex,btManifoldPoint& cp,
                                                            btScalar combinedTorsionalFriction,
                                                            btCollisionObject* colObj0,btCollisionObject* colObj1, btScalar relaxation, const btContactSolverInfo& infoGlobal, btScalar desiredVelocity=0, btScalar cfmSlip=0);

	void setupMultiBodyJointLimitConstraint(btMultiBodySolverConstraint& constraintRow, 
																 btScalar* jacA,btScalar* jacB,
																 btScalar penetration,btScalar combinedFrictionCoeff, btScalar combinedRestitutionCoeff,
																 const btContactSolverInfo& infoGlobal);

	void setupMultiBodyContactConstraint(btMultiBodySolverConstraint& solverConstraint, 
																 const btVector3& contactNormal,
																 btManifoldPoint& cp, const btContactSolverInfo& infoGlobal,
																 btScalar& relaxation,
																 bool isFriction, btScalar desiredVelocity=0, btScalar cfmSlip=0);
    
    //either rolling or spinning friction
    void setupMultiBodyTorsionalFrictionConstraint(btMultiBodySolverConstraint& solverConstraint,
                                         const btVector3& contactNormal,
                                         btManifoldPoint& cp,
                                        btScalar combinedTorsionalFriction,
                                        const btContactSolverInfo& infoGlobal,
                                         btScalar& relaxation,
                                         bool isFriction, btScalar desiredVelocity=0, btScalar cfmSlip=0);

	void convertMultiBodyContact(btPersistentManifold* manifold,const btContactSolverInfo& infoGlobal);
	virtual btScalar solveGroupCacheFriendlySetup(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
//	virtual btScalar solveGroupCacheFriendlyIterations(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);

	virtual btScalar solveSingleIteration(int iteration, btCollisionObject** bodies ,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
	void	applyDeltaVee(btScalar* deltaV, btScalar impulse, int velocityIndex, int ndof);
	void writeBackSolverBodyToMultiBody(btMultiBodySolverConstraint& constraint, btScalar deltaTime);
public:

	struct btMultiBodyConstraints
	{
		btSequentialImpulseConstraintSolver::btConstraints m_rigidBodyConstraints;

		/// Multibody (joint) constraints. This is shared by all the blocks.
		btMultiBodyConstraint** m_multiBodyConstraints;

		/// Number of multibody (joint) constraints. This is shared by all the
		/// blocks.
		int m_numMultiBodyConstraints;

		/// Array of multibody non-contact constraints
		btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyNonContactConstraints;

		/// Array of multibody normal contact constraints
		btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyNormalContactConstraints;

		/// Array of multibody friction contact constraints
		btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyFrictionContactConstraints;

		/// Array of multibody rolling friction contact constraints
		btAlignedObjectArray<btMultiBodySolverConstraint> m_multiBodyTorsionalFrictionContactConstraints;

		/// Pointer to the block constraint solver's multi body Jacobian data, which
		/// is shared by all the constraint blocks.
		btMultiBodyJacobianData* m_data;
	};

	BT_DECLARE_ALIGNED_ALLOCATOR();

	virtual void setMultiBodyConstraints(const btMultiBodyConstraints& data);

	///this method should not be called, it was just used during porting/integration of Featherstone btMultiBody, providing backwards compatibility but no support for btMultiBodyConstraint (only contact constraints)
	virtual btScalar solveGroup(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifold,int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& info, btIDebugDraw* debugDrawer,btDispatcher* dispatcher);
	virtual btScalar solveGroupConvertConstraintPrestep(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
	virtual btScalar solveGroupConvertConstraintPoststep(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
	virtual btScalar solveSingleIterationNew(int iteration, btCollisionObject** bodies ,int numBodies,btPersistentManifold** manifoldPtr, int numManifolds,btTypedConstraint** constraints,int numConstraints,const btContactSolverInfo& infoGlobal,btIDebugDraw* debugDrawer);
	virtual btScalar solveGroupCacheFriendlyFinish(btCollisionObject** bodies,int numBodies,const btContactSolverInfo& infoGlobal);
	
	virtual void solveMultiBodyGroup(btCollisionObject** bodies,int numBodies,btPersistentManifold** manifold,int numManifolds,btTypedConstraint** constraints,int numConstraints,btMultiBodyConstraint** multiBodyConstraints, int numMultiBodyConstraints, const btContactSolverInfo& info, btIDebugDraw* debugDrawer,btDispatcher* dispatcher);
};

	
	


#endif //BT_MULTIBODY_CONSTRAINT_SOLVER_H

