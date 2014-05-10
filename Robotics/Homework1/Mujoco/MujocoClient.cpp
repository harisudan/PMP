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
void setGrip(double gripAmount, Vector &thetahat)
{
  // indices for joints of two gripper fingers
  static const int L_GRIP_JOINT_INDEX = 7;
  static const int R_GRIP_JOINT_INDEX = 8;
  
  thetahat[L_GRIP_JOINT_INDEX] = -gripAmount;
  thetahat[R_GRIP_JOINT_INDEX] = +gripAmount;
}

void main(void)
{
  // indices for nodes in the scene
	static const int TARGET_GEOM_INDEX = 1;
	static const int HAND_GEOM_INDEX = 4;
	static const int OBJECT_GEOM_INDEX = 13; 

	// connect to mujoco server
	mjInit();
	mjConnect(10000);

	// load hand model
	if(mjIsConnected())
	{
		mjLoad("hand.xml");
		Sleep(1000);  // wait till the load is complete
	}

	if(mjIsConnected() && mjIsModel())
	{
		mjSetMode(2);
		mjReset();

		// size containts model dimensions
		mjSize size = mjGetSize();

		// number of arm degress of freedom (DOFs) and controls
		// (does not include target object degrees of freedom)
		int dimtheta = size.nu;

		// identity matrix (useful for implementing Jacobian pseudoinverse methods)
		Matrix I(dimtheta, dimtheta);
		I.setIdentity();

		// identity matrix (useful for implementing Jacobian pseudoinverse methods)
		Matrix I3By3(3, 3);
		I3By3.setIdentity();

		// Zero Reference vector
		Vector ZeroVector(dimtheta);
		ZeroVector.setConstant(0.0);

		// target arm DOFs
		Vector thetahat(dimtheta);
		thetahat.setConstant(0.0);

		for(;;) // run simulation forever
		{
			// simulation advance substep
			mjStep1();

			mjState state = mjGetState();
			// current arm degrees of freedom
			Vector theta(dimtheta);
			// state.qpos contains DOFs for the whole system. Here extract
			// only DOFs for the arm into theta			
			for(int e=0; e<dimtheta; e++)
			{
				theta[e] = state.qpos[e];
			}

			// Original provided code- that was used for all HW questions except 4
			// 
			//mjCartesian geoms = mjGetGeoms();
			//// current hand position
			//Vector x(3, geoms.pos + 3*HAND_GEOM_INDEX);
			//// target hand position
			//Vector xhat(3, geoms.pos + 3*TARGET_GEOM_INDEX);
			//// current hand orientation
			//Quat r(geoms.quat + 4*HAND_GEOM_INDEX);
			//// target hand orientation
			//Quat rhat(geoms.quat + 4*TARGET_GEOM_INDEX);

			//mjJacobian jacobians = mjJacGeom(HAND_GEOM_INDEX);
			//// current hand position Jacobian
			//Matrix Jpos(3,jacobians.nv, jacobians.jacpos);
			//// current hand orientation Jacobian
			//Matrix Jrot(3,jacobians.nv, jacobians.jacrot);
			//// extract only columns of the Jacobian that relate to arm DOFs
			//Jpos = Jpos.getBlock(0,3, 0, dimtheta);
			//Jrot = Jrot.getBlock(0,3, 0, dimtheta);

			mjCartesian geoms = mjGetGeoms();
			// current hand position
			Vector x(3, geoms.pos + 3*HAND_GEOM_INDEX);
			// target object position - which has the object to pickup
			Vector xhatObject(3, geoms.pos + 3*OBJECT_GEOM_INDEX);
			// target hand position - which has the target position
			Vector xhat(3, geoms.pos + 3*TARGET_GEOM_INDEX);
			// current hand orientation
			Quat r(geoms.quat + 4*HAND_GEOM_INDEX);
			// target object orientation
			// Quat rhatObject(geoms.quat + 4 * OBJECT_GEOM_INDEX);
			// target hand orientation
			Quat rhat(geoms.quat + 4*TARGET_GEOM_INDEX);

			mjJacobian jacobians = mjJacGeom(HAND_GEOM_INDEX);
			// current hand position Jacobian
			Matrix Jpos(3,jacobians.nv, jacobians.jacpos);
			// current hand orientation Jacobian
			Matrix Jrot(3,jacobians.nv, jacobians.jacrot);
			// extract only columns of the Jacobian that relate to arm DOFs
			Jpos = Jpos.getBlock(0,3, 0, dimtheta);
			Jrot = Jrot.getBlock(0,3, 0, dimtheta);

			// -- your code goes here --
			// set thetahat through Jacobian control methods here
			// for part 4 of assignment, you may need to call setGrip(amount, thetahat) here
			// thetahat should be filled by this point

			// Jacobian Transponse Method.
			// 
			/*Vector xdelta = xhat - x;
			Matrix Jpostranspose = Jpos.transpose();
			float alpha = 0.01;
			Vector deltatheta = (Jpostranspose.operator*(xdelta)) * alpha;
			thetahat = theta + deltatheta;*/

			// Pseudo Inverse Method - Basic
			// 
			//Vector xdelta = xhat - x;
			//Matrix Jpostransponse = Jpos.transpose();
			//float alpha = 0.01;
			//Matrix JJpostransponseInverse = (Jpos * Jpostransponse).inverse();
			//Vector deltaTheta = Jpostransponse * (JJpostransponseInverse * xdelta);	// Multiplied with xdelta to perform Matrix * vector rather than Matrix * Matrix 
			//deltaTheta = deltaTheta * alpha;
			//thetahat = deltaTheta + theta;

			// Pseudo Inverse Method - Explicit Optimization Criterion
			// 
			/*Vector xdelta = xhat - x;
			Matrix Jpostransponse = Jpos.transpose();
			float alpha = 0.1;
			Matrix JHash = Jpostransponse * ((Jpos * Jpostransponse).inverse());
			Matrix JHashJ = JHash * Jpos;
			Matrix Intermediate1 = I - JHashJ;
			Vector Intermediate2 = ZeroVector - theta;
			Vector RightTerm = Intermediate1 * Intermediate2;
			Vector LeftTerm = (JHash * xdelta) * alpha;
			Vector deltaTheta = LeftTerm + RightTerm;
			thetahat = deltaTheta + theta;*/

			// Jacobian Transponse for orientation control
			// 
			/*Vector rDelta = quatdiff(r, rhat);
			Matrix Jrottransponse = Jrot.transpose();
			float alpha2 = 0.01;
			Vector deltaTheta = (Jrottransponse * rDelta) * alpha2;
			thetahat = theta + deltaTheta;*/

			// Pseudo Inverse Method for orientation control - Basic
			// 
			//Vector rDelta = quatdiff(r, rhat);
			//Matrix Jrottransponse = Jrot.transpose();
			//float alpha2 = 0.01;
			//Matrix JJrottransponseInverse = (Jrot * Jrottransponse).inverse();
			//Vector deltaTheta = Jrottransponse * (JJrottransponseInverse * rDelta);	// Multiplied with xdelta to perform Matrix * vector rather than Matrix * Matrix 
			//deltaTheta = deltaTheta * alpha2;
			//thetahat = deltaTheta + theta;

			// Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
			// 
			//Vector rdelta = quatdiff(r, rhat);
			//Matrix Jrottransponse = Jrot.transpose();
			//float alpha2 = 0.1;
			//Matrix JHash = Jrottransponse * ((Jrot * Jrottransponse).inverse());
			//Matrix JHashJ = JHash * Jrot;
			//Matrix Intermediate1 = I - JHashJ;
			//Vector Intermediate2 = ZeroVector - theta;
			//Vector RightTerm = Intermediate1 * Intermediate2;
			//Vector LeftTerm = (JHash * rdelta) * alpha2;
			//Vector deltaTheta = LeftTerm + RightTerm;
			//thetahat = deltaTheta + theta; 

			// Common code for concatenating position and orientation variables
			// Copy used for HW Part - 3
			// 
			/*Vector rDelta = quatdiff(r, rhat);
			Vector xdelta = xhat - x;
			Vector ConcatenatedXRDelta(xdelta.getSize() + rDelta.getSize());
			ConcatenatedXRDelta[0] = xdelta[0];
			ConcatenatedXRDelta[1] = xdelta[1];
			ConcatenatedXRDelta[2] = xdelta[2];
			ConcatenatedXRDelta[3] = rDelta[0];
			ConcatenatedXRDelta[4] = rDelta[1];
			ConcatenatedXRDelta[5] = rDelta[2];
			assert(Jrot.getNumCols() == Jpos.getNumCols());
			assert(Jrot.getNumRows() == Jpos.getNumRows());
			Matrix ConcatenatedJacobMatrix(Jpos.getNumRows() + Jrot.getNumRows(), Jpos.getNumCols());
			
			int i = 0, j = 0;
			for (i = 0; i < Jpos.getNumRows(); i++)
			{
				ConcatenatedJacobMatrix.setRow(i, Jpos.getRow(i));
			}

			assert(i == Jpos.getNumRows());

			for (j = 0; j < Jrot.getNumRows(); j++)
			{
				ConcatenatedJacobMatrix.setRow(i + j, Jrot.getRow(j));
			}

			assert(j == Jrot.getNumRows());*/

			// Position and orientation - Jacobian Transpose method
			// 
			/*Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
			float alpha3 = 0.01;
			Vector deltatheta = (JConcatTranspose * ConcatenatedXRDelta) * alpha3;
			thetahat = theta + deltatheta;*/

			// Position and orientation - PseudoInverse Method - basic
			//
			//Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
			//float alpha4 = 0.01;
			//Matrix JJConcatTransponseInverse = (ConcatenatedJacobMatrix * JConcatTranspose).inverse();
			//Vector deltaTheta = JConcatTranspose * (JJConcatTransponseInverse * ConcatenatedXRDelta);	// Multiplied with xdelta to perform Matrix * vector rather than Matrix * Matrix 
			//deltaTheta = deltaTheta * alpha4;
			//thetahat = deltaTheta + theta;

			// Position and orientation - PseudoInverse with Explicit Optimization Criterion
			//
			/*Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
			float alpha5 = 0.01;
			Matrix JHash = JConcatTranspose * ((ConcatenatedJacobMatrix * JConcatTranspose).inverse());
			Matrix JHashJ = JHash * ConcatenatedJacobMatrix;
			Matrix Intermediate1 = I - JHashJ;
			Vector Intermediate2 = ZeroVector - theta;
			Vector RightTerm = Intermediate1 * Intermediate2;
			Vector LeftTerm = (JHash * ConcatenatedXRDelta) * alpha5;
			Vector deltaTheta = LeftTerm + RightTerm;
			thetahat = deltaTheta + theta;*/

			// set target DOFs to thetahat and advance simulation
			mjSetControl(dimtheta, thetahat);
			mjStep2();
		}

		mjDisconnect();
	}
	mjClear();
}
