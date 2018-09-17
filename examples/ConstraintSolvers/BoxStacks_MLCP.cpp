#include "BoxStacks_MLCP.h"
#include "../OpenGLWindow/SimpleOpenGL3App.h"
#include "btBulletDynamicsCommon.h"

#include "BulletDynamics/MLCPSolvers/btDantzigSolver.h"
#include "BulletDynamics/MLCPSolvers/btSolveProjectedGaussSeidel.h"

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

#include "../OpenGLWindow/GLInstancingRenderer.h"
#include "BulletCollision/CollisionShapes/btShapeHull.h"

#include "../CommonInterfaces/CommonMultiBodyBase.h"

#include "../RobotSimulator/b3RobotSimulatorClientAPI.h"
#include "../Importers/ImportURDFDemo/BulletUrdfImporter.h"
#include "../Importers/ImportURDFDemo/MyMultiBodyCreator.h"
#include "../Importers/ImportURDFDemo/URDF2Bullet.h"

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

#define ARRAY_SIZE_Y 4
#define ARRAY_SIZE_X 4
#define ARRAY_SIZE_Z 4

class BoxStacks_MLCP : public CommonMultiBodyBase
{
public:
	BoxStacks_MLCP(GUIHelperInterface* helper);
	virtual ~BoxStacks_MLCP();

	virtual void initPhysics();

	virtual void stepSimulation(float deltaTime);

	virtual void resetCamera()
	{
		float dist = 1;
		float pitch = -35;
		float yaw = 50;
		float targetPos[3] = {-3, 2.8, -2.5};
		m_guiHelper->resetCamera(dist, yaw, pitch, targetPos[0], targetPos[1], targetPos[2]);
	}

	void createBoxStack(int numBoxes, btScalar centerX, btScalar centerY);
	btMultiBody* createFeatherstoneMultiBody(class btMultiBodyDynamicsWorld* world, int numLinks, const btVector3& basePosition, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents, bool spherical = false, bool fixedBase = false);
	void createGround(const btVector3& halfExtents = btVector3(50, 50, 50), btScalar zOffSet = btScalar(-1.55));
	void addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents);
	btMultiBody* loadRobot(std::string filepath = "kuka_iiwa/model.urdf");
};

static bool g_fixedBase = true;
static bool g_firstInit = true;
static float scaling = 0.4f;
static float friction = 1.;
static int g_constraintSolverType = 0;

BoxStacks_MLCP::BoxStacks_MLCP(GUIHelperInterface* helper)
	: CommonMultiBodyBase(helper)
{
	m_guiHelper->setUpAxis(1);
}

BoxStacks_MLCP::~BoxStacks_MLCP()
{
	// Do nothing
}

void BoxStacks_MLCP::stepSimulation(float deltaTime)
{
	//use a smaller internal timestep, there are stability issues
	float internalTimeStep = 1. / 240.f;
	m_dynamicsWorld->stepSimulation(deltaTime, 10, internalTimeStep);
}

void BoxStacks_MLCP::initPhysics()
{
	m_guiHelper->setUpAxis(1);

	createEmptyDynamicsWorld();




	if (g_constraintSolverType == 3)
	{
		g_constraintSolverType = 0;
	}

	btMultiBodyConstraintSolver* sol = new btMultiBodyBlockConstraintSolver();
	btMLCPSolverInterface* mlcp;
	switch (g_constraintSolverType++)
	{
		case 0:
			sol = new btMultiBodyConstraintSolver;
			b3Printf("Constraint Solver: Sequential Impulse");
			break;
		case 1:
			mlcp = new btSolveProjectedGaussSeidel();
			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
			b3Printf("Constraint Solver: MLCP + PGS");
			break;
		case 2:
			mlcp = new btDantzigSolver();
			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
			b3Printf("Constraint Solver: MLCP + Dantzig");
			break;
//		default:
//			mlcp = new btLemkeSolver();
//			sol = new btMultiBodyMLCPConstraintSolver(mlcp);
//			b3Printf("Constraint Solver: MLCP + Lemke");
//			break;
	}
	m_solver = sol;

	btMultiBodyDynamicsWorld* world = new btMultiBodyDynamicsWorld(m_dispatcher, m_broadphase, sol, m_collisionConfiguration);
	m_dynamicsWorld = world;





	//m_dynamicsWorld->setGravity(btVector3(0,0,0));
	m_guiHelper->createPhysicsDebugDrawer(m_dynamicsWorld);

	if (m_dynamicsWorld->getDebugDrawer())
		m_dynamicsWorld->getDebugDrawer()->setDebugMode(btIDebugDraw::DBG_DrawWireframe+btIDebugDraw::DBG_DrawContactPoints);

	///create a few basic rigid bodies
	btBoxShape* groundShape = createBoxShape(btVector3(btScalar(50.),btScalar(50.),btScalar(50.)));


	//groundShape->initializePolyhedralFeatures();
	//btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0,1,0),50);

	m_collisionShapes.push_back(groundShape);

	btTransform groundTransform;
	groundTransform.setIdentity();
	groundTransform.setOrigin(btVector3(0,-50,0));

	{
		btScalar mass(0.);
		createRigidBody(mass,groundTransform,groundShape, btVector4(0,0,1,1));
	}


	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btBoxShape* colShape = createBoxShape(btVector3(.1,.1,.1));


		//btCollisionShape* colShape = new btSphereShape(btScalar(1.));
		m_collisionShapes.push_back(colShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		btScalar	mass(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
			colShape->calculateLocalInertia(mass,localInertia);


		for (int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for (int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					startTransform.setOrigin(btVector3(
										btScalar(0.2*i),
										btScalar(2+.2*k),
										btScalar(0.2*j)));


					createRigidBody(mass,startTransform,colShape);


				}
			}
		}
	}


	m_guiHelper->autogenerateGraphicsObjects(m_dynamicsWorld);
}

void BoxStacks_MLCP::createBoxStack(int numBoxes, btScalar centerX, btScalar centerZ)
{
	//create a few dynamic rigidbodies
	// Re-using the same collision is better for memory usage and performance

	const btScalar boxHalfSize = btScalar(0.1);

	btBoxShape* colShape = createBoxShape(btVector3(boxHalfSize, boxHalfSize, boxHalfSize));
	m_collisionShapes.push_back(colShape);

	/// Create Dynamic Objects
	btTransform startTransform;
	startTransform.setIdentity();

	btScalar mass(1.0);

	btVector3 localInertia(0, 0, 0);
	colShape->calculateLocalInertia(mass,localInertia);

	for (int i = 0; i < numBoxes; ++i)
	{
		startTransform.setOrigin(btVector3(centerX, btScalar(btScalar(2) * boxHalfSize * i), centerZ));
		createRigidBody(mass, startTransform, colShape);
	}
}

btMultiBody* BoxStacks_MLCP::createFeatherstoneMultiBody(btMultiBodyDynamicsWorld* pWorld, int numLinks, const btVector3& basePosition, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents, bool spherical, bool fixedBase)
{
	//init the base
	btVector3 baseInertiaDiag(0.f, 0.f, 0.f);
	float baseMass = 1.f;

	if (baseMass)
	{
		btCollisionShape* pTempBox = new btBoxShape(btVector3(baseHalfExtents[0], baseHalfExtents[1], baseHalfExtents[2]));
		pTempBox->calculateLocalInertia(baseMass, baseInertiaDiag);
		delete pTempBox;
	}

	bool canSleep = false;

	btMultiBody* pMultiBody = new btMultiBody(numLinks, baseMass, baseInertiaDiag, fixedBase, canSleep);

	btQuaternion baseOriQuat(0.f, 0.f, 0.f, 1.f);
	pMultiBody->setBasePos(basePosition);
	pMultiBody->setWorldToBaseRot(baseOriQuat);
	btVector3 vel(0, 0, 0);

	//init the links
	btVector3 hingeJointAxis(1, 0, 0);
	float linkMass = 1.f;
	btVector3 linkInertiaDiag(0.f, 0.f, 0.f);

	btCollisionShape* pTempBox = new btBoxShape(btVector3(linkHalfExtents[0], linkHalfExtents[1], linkHalfExtents[2]));
	pTempBox->calculateLocalInertia(linkMass, linkInertiaDiag);
	delete pTempBox;

	//y-axis assumed up
	btVector3 parentComToCurrentCom(0, -linkHalfExtents[1] * 2.f, 0);                      //par body's COM to cur body's COM offset
	btVector3 currentPivotToCurrentCom(0, -linkHalfExtents[1], 0);                         //cur body's COM to cur body's PIV offset
	btVector3 parentComToCurrentPivot = parentComToCurrentCom - currentPivotToCurrentCom;  //par body's COM to cur body's PIV offset

	//////
	btScalar q0 = 0.f * SIMD_PI / 180.f;
	btQuaternion quat0(btVector3(0, 1, 0).normalized(), q0);
	quat0.normalize();
	/////

	for (int i = 0; i < numLinks; ++i)
	{
		if (!spherical)
			pMultiBody->setupRevolute(i, linkMass, linkInertiaDiag, i - 1, btQuaternion(0.f, 0.f, 0.f, 1.f), hingeJointAxis, parentComToCurrentPivot, currentPivotToCurrentCom, true);
		else
			//pMultiBody->setupPlanar(i, linkMass, linkInertiaDiag, i - 1, btQuaternion(0.f, 0.f, 0.f, 1.f)/*quat0*/, btVector3(1, 0, 0), parentComToCurrentPivot*2, false);
			pMultiBody->setupSpherical(i, linkMass, linkInertiaDiag, i - 1, btQuaternion(0.f, 0.f, 0.f, 1.f), parentComToCurrentPivot, currentPivotToCurrentCom, true);
	}

	pMultiBody->finalizeMultiDof();

	///
	pWorld->addMultiBody(pMultiBody);
	///
	return pMultiBody;
}

void BoxStacks_MLCP::createGround(const btVector3& halfExtents, btScalar zOffSet)
{
	btCollisionShape* groundShape = new btBoxShape(halfExtents);
	m_collisionShapes.push_back(groundShape);

	// rigidbody is dynamic if and only if mass is non zero, otherwise static
	btScalar mass(0.);
	const bool isDynamic = (mass != 0.f);

	btVector3 localInertia(0, 0, 0);
	if (isDynamic)
		groundShape->calculateLocalInertia(mass, localInertia);

	// using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
	btTransform groundTransform;
	groundTransform.setIdentity();
	groundTransform.setOrigin(btVector3(0, -halfExtents.z() + zOffSet, 0));
	btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform);
	btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, groundShape, localInertia);
	btRigidBody* body = new btRigidBody(rbInfo);

	// add the body to the dynamics world
	m_dynamicsWorld->addRigidBody(body, 1, 1 + 2);
}

void BoxStacks_MLCP::addColliders(btMultiBody* pMultiBody, btMultiBodyDynamicsWorld* pWorld, const btVector3& baseHalfExtents, const btVector3& linkHalfExtents)
{
	btAlignedObjectArray<btQuaternion> world_to_local;
	world_to_local.resize(pMultiBody->getNumLinks() + 1);

	btAlignedObjectArray<btVector3> local_origin;
	local_origin.resize(pMultiBody->getNumLinks() + 1);
	world_to_local[0] = pMultiBody->getWorldToBaseRot();
	local_origin[0] = pMultiBody->getBasePos();

	{
		btScalar quat[4] = {-world_to_local[0].x(), -world_to_local[0].y(), -world_to_local[0].z(), world_to_local[0].w()};

		if (1)
		{
			btCollisionShape* box = new btBoxShape(baseHalfExtents);
			btMultiBodyLinkCollider* col = new btMultiBodyLinkCollider(pMultiBody, -1);
			col->setCollisionShape(box);

			btTransform tr;
			tr.setIdentity();
			tr.setOrigin(local_origin[0]);
			tr.setRotation(btQuaternion(quat[0], quat[1], quat[2], quat[3]));
			col->setWorldTransform(tr);

			pWorld->addCollisionObject(col, 2, 1 + 2);

			col->setFriction(friction);
			pMultiBody->setBaseCollider(col);
		}
	}

	for (int i = 0; i < pMultiBody->getNumLinks(); ++i)
	{
		const int parent = pMultiBody->getParent(i);
		world_to_local[i + 1] = pMultiBody->getParentToLocalRot(i) * world_to_local[parent + 1];
		local_origin[i + 1] = local_origin[parent + 1] + (quatRotate(world_to_local[i + 1].inverse(), pMultiBody->getRVector(i)));
	}

	for (int i = 0; i < pMultiBody->getNumLinks(); ++i)
	{
		btVector3 posr = local_origin[i + 1];

		btScalar quat[4] = {-world_to_local[i + 1].x(), -world_to_local[i + 1].y(), -world_to_local[i + 1].z(), world_to_local[i + 1].w()};

		btCollisionShape* box = new btBoxShape(linkHalfExtents);
		btMultiBodyLinkCollider* col = new btMultiBodyLinkCollider(pMultiBody, i);

		col->setCollisionShape(box);
		btTransform tr;
		tr.setIdentity();
		tr.setOrigin(posr);
		tr.setRotation(btQuaternion(quat[0], quat[1], quat[2], quat[3]));
		col->setWorldTransform(tr);
		col->setFriction(friction);
		pWorld->addCollisionObject(col, 2, 1 + 2);

		pMultiBody->getLink(i).m_collider = col;
	}
}

btMultiBody*BoxStacks_MLCP::loadRobot(std::string filepath)
{
	btMultiBody* m_multiBody = 0;
	BulletURDFImporter u2b(m_guiHelper,0,1,0);
	bool loadOk = u2b.loadURDF(filepath.c_str());// lwr / kuka.urdf");
	if (loadOk)
	{
			int rootLinkIndex = u2b.getRootLinkIndex();
			b3Printf("urdf root link index = %d\n",rootLinkIndex);
			MyMultiBodyCreator creation(m_guiHelper);
			btTransform identityTrans;
			identityTrans.setIdentity();
			ConvertURDF2Bullet(u2b,creation, identityTrans,m_dynamicsWorld,true,u2b.getPathPrefix());
			for (int i = 0; i < u2b.getNumAllocatedCollisionShapes(); i++)
			{
				m_collisionShapes.push_back(u2b.getAllocatedCollisionShape(i));
			}
			m_multiBody = creation.getBulletMultiBody();
			if (m_multiBody)
			{
				//kuka without joint control/constraints will gain energy explode soon due to timestep/integrator
				//temporarily set some extreme damping factors until we have some joint control or constraints
				m_multiBody->setAngularDamping(0*0.99);
				m_multiBody->setLinearDamping(0*0.99);
				b3Printf("Root link name = %s",u2b.getLinkName(u2b.getRootLinkIndex()).c_str());
			}
	}
	return m_multiBody;
}

CommonExampleInterface* BoxStacks_MLCPCreateFunc(CommonExampleOptions& options)
{
	return new BoxStacks_MLCP(options.m_guiHelper);
}
