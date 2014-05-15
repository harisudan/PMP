// mujoco communication interface
#include "comm.h"

// algebra classes
#include "GMatrix.h"
#include "GVector.h"
#include "Quaternion.h"
#include "windows.h"

typedef GVector<double> Vector;
typedef GMatrix<double> Matrix;
typedef Quaternion<double> Quat;

// sets thetahat such that gripper fingers open or close by gripAmount
// gripAmount = 0: default gipper posture
// gripAmount > 0: gripper fingers opened
// gripAmount < 0: gripper fingers closed
// between -1 and +1 are acceptable gripAmount values
void setGrip2(double gripAmount, Vector &thetahat)
{
	// indices for joints of two gripper fingers
	static const int L_GRIP_JOINT_INDEX = 7;
	static const int R_GRIP_JOINT_INDEX = 8;

	thetahat[L_GRIP_JOINT_INDEX] = -gripAmount;
	thetahat[R_GRIP_JOINT_INDEX] = +gripAmount;
}

void main2(void)
{
	// indices for nodes in the scene
	static const int TARGET_GEOM_INDEX = 1;
	static const int HAND_GEOM_INDEX = 4;
	static const int OBJECT_GEOM_INDEX = 13;

	// connect to mujoco server.
	// 
	mjInit();
	mjConnect(10000);

	double lTarget = 0.0;
	double rTarget = 0.0;

	// load hand model.
	// 
	if (mjIsConnected())
	{
		mjLoad("hand.xml");
		Sleep(1000);  // wait till the load is complete
	}

	if (mjIsConnected() && mjIsModel())
	{
		mjSetMode(2);
		mjReset();

		// size containts model dimensions.
		// 
		mjSize size = mjGetSize();

		// number of arm degress of freedom (DOFs) and controls
		// (does not include target object degrees of freedom)
		// 
		int dimtheta = size.nu;

		// identity matrix (useful for implementing Jacobian pseudoinverse methods).
		// 
		Matrix I(dimtheta, dimtheta);
		I.setIdentity();

		// Zero Reference vector.
		// 
		Vector ZeroVector(dimtheta);
		ZeroVector.setConstant(0.0);

		// target arm DOFs.
		// 
		Vector thetahat(dimtheta);		thetahat.setConstant(0.0);

		for (;;) // run simulation forever
		{
			// simulation advance substep.
			// 
			mjStep1();

			mjState state = mjGetState();

			// current arm degrees of freedom.
			// 
			Vector theta(dimtheta);

			// state.qpos contains DOFs for the whole system. Here extract
			// only DOFs for the arm into theta.
			// 
			for (int e = 0; e<dimtheta; e++)
			{
				theta[e] = state.qpos[e];
			}

			mjCartesian geoms = mjGetGeoms();

			// current hand/gripper position.
			// 
			Vector x(3, geoms.pos + 3 * HAND_GEOM_INDEX);

			// Target blue object position - which has the object to pickup
			// 
			Vector xhatObject(3, geoms.pos + 3 * OBJECT_GEOM_INDEX);

			// Target Green marker position.
			// 
			Vector xhat(3, geoms.pos + 3 * TARGET_GEOM_INDEX);

			// current hand/gripper orientation
			// 
			Quat r(geoms.quat + 4 * HAND_GEOM_INDEX);

			// Blue Object orientation
			// 
			Quat rhatObject(geoms.quat + 4 * OBJECT_GEOM_INDEX);

			// Target Green Marker orientation.
			// 
			Quat rhat(geoms.quat + 4 * TARGET_GEOM_INDEX);

			mjJacobian jacobians = mjJacGeom(HAND_GEOM_INDEX);

			// current hand position Jacobian.
			// 
			Matrix Jpos(3, jacobians.nv, jacobians.jacpos);

			// current hand orientation Jacobian.
			// 
			Matrix Jrot(3, jacobians.nv, jacobians.jacrot);

			// extract only columns of the Jacobian that relate to arm DOFs
			// 
			Jpos = Jpos.getBlock(0, 3, 0, dimtheta);
			Jrot = Jrot.getBlock(0, 3, 0, dimtheta);

			// -- your code goes here --
			// set thetahat through Jacobian control methods here
			// for part 4 of assignment, you may need to call setGrip(amount, thetahat) here
			// thetahat should be filled by this point			

			// ----------------------------------------------------------------------------------------
			// Homework Part 2
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// Orientation Control - Jacobian Transponse Method.
			// ----------------------------------------------------------------------------------------
			/*Vector rDelta = quatdiff(r, rhat);
			Matrix Jrottransponse = Jrot.transpose();
			float alpha2 = 0.0001;
			Vector deltaTheta = (Jrottransponse * rDelta) * alpha2;
			thetahat = theta + deltaTheta;*/

			// ----------------------------------------------------------------------------------------
			// Orientation Control - Pseudo Inverse Method - Basic Method (Lecture Slide Page 10)
			// ----------------------------------------------------------------------------------------
			/*Vector rDelta = quatdiff(r, rhat);
			Matrix Jrottransponse = Jrot.transpose();
			float alpha2 = 0.01;
			Matrix JJrottransponseInverse = (Jrot * Jrottransponse).inverse();
			Vector deltaTheta = Jrottransponse * (JJrottransponseInverse * rDelta);	// Multiplied with xdelta to perform Matrix * vector rather than Matrix * Matrix
			deltaTheta = deltaTheta * alpha2;
			thetahat = deltaTheta + theta;*/

			// ----------------------------------------------------------------------------------------
			// Orientation Control - Pseudo Inverse Method - Explicit Optimization Criterion (Lecture Slide Page 13)
			// ----------------------------------------------------------------------------------------
			Vector rdelta = quatdiff(r, rhat);
			Matrix Jrottransponse = Jrot.transpose();
			float alpha2 = 0.01;
			Matrix JHash = Jrottransponse * ((Jrot * Jrottransponse).inverse());
			Matrix JHashJ = JHash * Jrot;
			Matrix Intermediate1 = I - JHashJ;
			Vector Intermediate2 = ZeroVector - theta;
			Vector RightTerm = Intermediate1 * Intermediate2;
			Vector LeftTerm = (JHash * rdelta) * alpha2;
			Vector deltaTheta = LeftTerm + RightTerm;
			thetahat = deltaTheta + theta;

			// set target DOFs to thetahat and advance simulation.
			// 
			mjSetControl(dimtheta, thetahat);
			mjStep2();
		}

		mjDisconnect();
	}
	mjClear();
}