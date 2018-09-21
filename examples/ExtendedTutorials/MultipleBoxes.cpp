/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2015 Google Inc. http://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/



#include "MultipleBoxes.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h" 
#include "../CommonInterfaces/CommonRigidBodyBase.h"
#include "../CommonInterfaces/CommonMultiBodyBase.h"

#include "BulletDynamics/Featherstone/btMultiBodyDynamicsWorld.h"
#include "BulletDynamics/Featherstone/btMultiBody.h"
#include "BulletDynamics/Featherstone/btMultiBodyConstraintSolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyBlockConstraintSolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyMLCPConstraintSolver.h"
#include "BulletDynamics/Featherstone/btMultiBodyDynamicsWorld.h"
#include "BulletDynamics/Featherstone/btMultiBodyLinkCollider.h"
#include "BulletDynamics/Featherstone/btMultiBodyLink.h"
#include "BulletDynamics/Featherstone/btMultiBodyJointLimitConstraint.h"
#include "BulletDynamics/Featherstone/btMultiBodyJointMotor.h"
#include "BulletDynamics/Featherstone/btMultiBodyPoint2Point.h"
#include "BulletDynamics/Featherstone/btMultiBodyFixedConstraint.h"
#include "BulletDynamics/Featherstone/btMultiBodySliderConstraint.h"

#include "BulletDynamics/MLCPSolvers/btDantzigSolver.h"
#include "BulletDynamics/MLCPSolvers/btLemkeSolver.h"
#include "BulletDynamics/MLCPSolvers/btSolveProjectedGaussSeidel.h"

static int g_constraintSolverType = 0;
const int TOTAL_BOXES = 10;
static int g_demo = 0;

struct MultipleBoxesExample : public CommonMultiBodyBase
{
	MultipleBoxesExample(struct GUIHelperInterface* helper)
		:CommonMultiBodyBase(helper)
	{
	}
	virtual ~MultipleBoxesExample(){}
	virtual void initPhysics();
	virtual void renderScene();
	void resetCamera()
	{
		float dist = 41;
		float pitch = -35;
		float yaw = 52;
		float targetPos[3]={0,0.46,0};
		m_guiHelper->resetCamera(dist,yaw,pitch,targetPos[0],targetPos[1],targetPos[2]);
	}
    void makeScene1();
    void makeScene2(btMultiBodyDynamicsWorld* world);
    btMultiBody* createSingleMultiBody(
            btMultiBodyDynamicsWorld* world,const char* name, float mass, const btTransform& startTransform, btCollisionShape* shape,  const btVector4& color = btVector4(1, 0, 0, 1));

};

btMultiBody* MultipleBoxesExample::createSingleMultiBody(
      btMultiBodyDynamicsWorld* world,
        const char* name, float mass, const btTransform& startTransform, btCollisionShape* shape,  const btVector4& color)
{
    btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));

    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);

    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        shape->calculateLocalInertia(mass, localInertia);

    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects

#define USE_MOTIONSTATE 1
#ifdef USE_MOTIONSTATE
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

    btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);

    btMultiBody* pMultiBody = new btMultiBody(0, mass, localInertia, false, false);
//    btTransform startTrans;
//    startTrans.setIdentity();
//    startTrans.setOrigin(btVector3(0,0,3));

    pMultiBody->setBaseWorldTransform(startTransform);

    btMultiBodyLinkCollider* col= new btMultiBodyLinkCollider(pMultiBody, -1);
    col->setCollisionShape(shape);
    pMultiBody->setBaseCollider(col);
    int collisionFilterGroup = isDynamic? int(btBroadphaseProxy::DefaultFilter) : int(btBroadphaseProxy::StaticFilter);
    int collisionFilterMask = isDynamic? 	int(btBroadphaseProxy::AllFilter) : 	int(btBroadphaseProxy::AllFilter ^ btBroadphaseProxy::StaticFilter);

    world->addCollisionObject(col,collisionFilterGroup,collisionFilterMask);//, 2,1+2);

    pMultiBody->finalizeMultiDof();
    pMultiBody->setBaseName(name);
    world->addMultiBody(pMultiBody);

    btAlignedObjectArray<btQuaternion> scratch_q;
    btAlignedObjectArray<btVector3> scratch_m;
    pMultiBody->forwardKinematics(scratch_q,scratch_m);
    btAlignedObjectArray<btQuaternion> world_to_local;
    btAlignedObjectArray<btVector3> local_origin;
    pMultiBody->updateCollisionObjectWorldTransforms(world_to_local,local_origin);

    //body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);

#else
    btRigidBody* body = new btRigidBody(mass, 0, shape, localInertia);
    body->setWorldTransform(startTransform);
#endif//

    return pMultiBody;
}

void MultipleBoxesExample::initPhysics()
{
	m_guiHelper->setUpAxis(1);

	createEmptyDynamicsWorld();
	
	if (g_constraintSolverType == 2)
	{
		g_constraintSolverType = 0;
	}

	g_demo = 1;

    if (g_demo == 0)
    {

        btMultiBodyConstraintSolver* sol = new btMultiBodyBlockConstraintSolver();
//        btMLCPSolverInterface* mlcp;
//        switch (g_constraintSolverType++)
//        {
//        case 0:
//            sol = new btMultiBodyConstraintSolver;
//            b3Printf("Constraint Solver: Sequential Impulse");
//            break;
//        case 1:
//            mlcp = new btSolveProjectedGaussSeidel();
//            sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//            b3Printf("Constraint Solver: MLCP + PGS");
//            break;
//        case 2:
//            mlcp = new btDantzigSolver();
//            sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//            b3Printf("Constraint Solver: MLCP + Dantzig");
//            break;
//        default:
//            mlcp = new btLemkeSolver();
//            sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//            b3Printf("Constraint Solver: MLCP + Lemke");
//            break;
//        }
		m_solver = sol;

		btMultiBodyDynamicsWorld* world = new btMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration);
		m_dynamicsWorld = world;

        auto& info = m_dynamicsWorld->getSolverInfo();
        info.m_globalCfm = 0.01;
		info.m_numIterations = 10;

    }

    btMultiBodyDynamicsWorld* world;
    if (g_demo == 1)
    {

        btMultiBodyConstraintSolver* sol = new btMultiBodyBlockConstraintSolver();
		btMLCPSolverInterface* mlcp;
		switch (g_constraintSolverType++)
		{
//		case 0:
//			sol = new btMultiBodyConstraintSolver;
//			b3Printf("Constraint Solver: Sequential Impulse");
//			break;
//		default:
//			sol = new btMultiBodyBlockConstraintSolver();
//			b3Printf("Constraint Solver: Block");
//			break;
//		case 1:
//			mlcp = new btSolveProjectedGaussSeidel();
//			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//			b3Printf("Constraint Solver: MLCP + PGS");
//			break;
//		default:
//			mlcp = new btDantzigSolver();
//			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//			b3Printf("Constraint Solver: MLCP + Dantzig");
//			break;
//		default:
//			mlcp = new btLemkeSolver();
//			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//			b3Printf("Constraint Solver: MLCP + Lemke");
//			break;
		}
        m_solver = sol;

        world = new btMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration);
        m_dynamicsWorld = world;

        auto& info = m_dynamicsWorld->getSolverInfo();
        info.m_globalCfm = 0.01;
		info.m_numIterations = 10;
    }

    m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);

	if (m_dynamicsWorld->getDebugDrawer())
		m_dynamicsWorld->getDebugDrawer()->setDebugMode(btIDebugDraw::DBG_DrawWireframe+btIDebugDraw::DBG_DrawContactPoints);

	///create a few basic rigid bodies
	btBoxShape* groundShape = createBoxShape(btVector3(btScalar(50.),btScalar(50.),btScalar(50.)));


//	m_collisionShapes.push_back(groundShape);

	btTransform groundTransform;
	groundTransform.setIdentity();
	groundTransform.setOrigin(btVector3(0,-50,0));
//	{
//		btScalar mass(0.);
//		createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1));
//	}

	createMultiBody("name", 0.0, groundTransform, groundShape, btVector4(0, 0, 1, 1));

    if (g_demo == 0)
        makeScene1();
    if (g_demo == 1)
        makeScene2(world);

	m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);
}


void MultipleBoxesExample::renderScene()
{
	CommonMultiBodyBase::renderScene();
}

void MultipleBoxesExample::makeScene1()
{

    //	btScalar maxMass = 10000;
    btScalar minMass = 0.1;
    btScalar stepMass = 1.5;

    {
        //create a few dynamic rigidbodies
        // Re-using the same collision is better for memory usage and performance
		//btBoxShape* colShape = createBoxShape(btVector3(1,1,1));
		btSphereShape* colShape = new btSphereShape(1);

        m_collisionShapes.push_back(colShape);

        /// Create Dynamic Objects
        btTransform startTransform;
        startTransform.setIdentity();

        btScalar	mass = minMass;

        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);

        btVector3 localInertia(0,0,0);
        if (isDynamic)
            colShape->calculateLocalInertia(mass,localInertia);


		//for(int i=0;i<TOTAL_BOXES;++i)
		for(int i=0;i<1;++i)
        {
            startTransform.setOrigin(btVector3(
                                         btScalar(0),
                                         btScalar(1+i*2),
                                         btScalar(0)));
            createRigidBody(mass,startTransform,colShape);
            mass *= stepMass;
        }
    }
}

void MultipleBoxesExample::makeScene2(btMultiBodyDynamicsWorld* world)
{
    btScalar distance = 5;

    //	btScalar maxMass = 10000;
    btScalar minMass = 0.1;
	btScalar stepMass = 2.0;

	const int half = TOTAL_BOXES/2;

    {
        //create a few dynamic rigidbodies
        // Re-using the same collision is better for memory usage and performance
		btBoxShape* colShape = createBoxShape(btVector3(1.5,1,1.5));

        m_collisionShapes.push_back(colShape);

        /// Create Dynamic Objects
        btTransform startTransform;
        startTransform.setIdentity();

        btScalar	mass = minMass;

        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);

        btVector3 localInertia(0,0,0);
        if (isDynamic)
            colShape->calculateLocalInertia(mass,localInertia);


		for(int i=0;i<half;++i)
		{
			startTransform.setOrigin(btVector3(
										 btScalar(-distance),
										 btScalar(2+i*2),
										 btScalar(0)));
			createSingleMultiBody(world, "group1a", mass,startTransform,colShape);
			mass *= stepMass;
		}
		for(int i=half;i<TOTAL_BOXES;++i)
		{
			startTransform.setOrigin(btVector3(
										 btScalar(-distance),
										 btScalar(2+i*2),
										 btScalar(0)));
			createSingleMultiBody(world, "group1b", mass,startTransform,colShape);
			mass *= stepMass;
		}
    }


    //	btScalar maxMass = 10000;
//    btScalar minMass = 0.1;
//    btScalar stepMass = 1.5;

    {
        //create a few dynamic rigidbodies
        // Re-using the same collision is better for memory usage and performance
		btBoxShape* colShape = createBoxShape(btVector3(1.5,1,1.5));

        m_collisionShapes.push_back(colShape);

        /// Create Dynamic Objects
        btTransform startTransform;
        startTransform.setIdentity();

        btScalar	mass = minMass;

        //rigidbody is dynamic if and only if mass is non zero, otherwise static
        bool isDynamic = (mass != 0.f);

        btVector3 localInertia(0,0,0);
        if (isDynamic)
            colShape->calculateLocalInertia(mass,localInertia);

		for(int i=0;i<half;++i)
		{
			startTransform.setOrigin(btVector3(
										 btScalar(distance),
										 btScalar(2+i*2),
										 btScalar(0)));
			createSingleMultiBody(world, "group2a", mass,startTransform,colShape);
			mass *= stepMass;
		}
		for(int i=half;i<TOTAL_BOXES;++i)
		{
			startTransform.setOrigin(btVector3(
										 btScalar(distance),
										 btScalar(2+i*2),
										 btScalar(0)));
			createSingleMultiBody(world, "group2b", mass,startTransform,colShape);
			mass *= stepMass;
		}
    }
}







CommonExampleInterface*    ET_MultipleBoxesCreateFunc(CommonExampleOptions& options)
{
	return new MultipleBoxesExample(options.m_guiHelper);
}



