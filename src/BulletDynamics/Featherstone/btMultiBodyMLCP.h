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

#ifndef BT_MULTIBODY_MLCP_H
#define BT_MULTIBODY_MLCP_H

#include "LinearMath/btMatrixX.h"
#include "LinearMath/btThreads.h"

class btMultiBodySolverConstraint;

/// Data struct for multibody MLCP block
struct btMultiBodyMLCP
{
	btMatrixXu m_A;
	btVectorXu m_b;
	//	btVectorXu m_bSplit;
	btVectorXu m_x;
	//	btVectorXu m_xSplit;
	btVectorXu m_lo;
	btVectorXu m_hi;
	btAlignedObjectArray<int> m_limitDependencies;
	btAlignedObjectArray<btMultiBodySolverConstraint*> m_allConstraintPtrArray;
};

#endif  // BT_MULTIBODY_MLCP_H
