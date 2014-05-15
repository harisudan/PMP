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

	enum StateMachine
	{
		Default,
		AlignPerpendicularToObject,
		ApproachingObject,
		OpenGrip,
		CloseGrip,
		ClosingGrip,
		MoveTowardsTarget,
	};

	StateMachine stateMachine = Default;

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

			// Create a quartenion which represents a rotation by 90 degrees around X-axis.
			// 
			Quat Rotation90(0.4, rhatObject[1] * 0.707, 0, 0);

			// Except when the state machine is in the last stage of taking the object towards the target marker,
			// mask the orientation of the BLUE object to be perpendicular to its actual orientation,
			// rotate around X-axis, by 90 degrees.
			// 
			if (stateMachine != MoveTowardsTarget)
			{
				rhatObject = Rotation90 * rhatObject;
			}

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
			// Homework Part 4
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// Step 1 - Align perpendicular to Capsule/Object orientation. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
			// ----------------------------------------------------------------------------------------
			// Step 2 - Move towards the BLUE object, with combined position and orientation control. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
			//			Orientation control remains perpendicular w.r.t BLUE object orientation.
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// Step 3 - When close enough to the BLUE object, open claws/set grip and continue combined position and orientation control. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
			//			Orientation control remains perpendicular w.r.t BLUE object orientation.
			// ----------------------------------------------------------------------------------------
			// Step 4 - Close grip, wait till closed to enough extent, based on the target close setting of -0.4
			// ----------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------
			// Step 5 - Take object, manouvre combined position and orientation control to move towards target GREEN marker.
			//			Orientation control such that BLUE object is oriented w.r.t GREEN marker.
			// ----------------------------------------------------------------------------------------

			Vector rdelta;
			rdelta = quatdiff(r, rhatObject);

			switch (stateMachine)
			{
				case Default:
					stateMachine = AlignPerpendicularToObject;
					break;
			}

			// If aligned perpendicularly is complete (which means rdelta[0] == 0), proceed to next state, ApproachingObject.
			// 
			if (rdelta[0] == 0.0 && stateMachine == AlignPerpendicularToObject)
			{
				stateMachine = ApproachingObject;
			}
			else if (stateMachine == AlignPerpendicularToObject)
			{
				// Continue with orientation control to complete step 1
				// 
				// ----------------------------------------------------------------------------------------
				// Step 1 - Align perpendicular to Capsule/Object orientation. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
				// ----------------------------------------------------------------------------------------
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
			}
			else if (stateMachine == ApproachingObject)
			{
				// ----------------------------------------------------------------------------------------
				// Step 2 - Move towards the BLUE object, with combined position and orientation control. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
				//			Orientation control remains perpendicular w.r.t BLUE object orientation.
				// ----------------------------------------------------------------------------------------
				
				Vector rDelta = quatdiff(r, rhatObject);
				Vector xdelta = xhatObject - x;
				
				// This threshold/constant is based on the size of the blue object
				// and the size of the claws/fingers
				// When close enough to the object, open the claws/gripper to grip the object.
				// 
				if (xdelta[0] < 0.11)
				{
					stateMachine = OpenGrip;
				}
				else
				{					
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

					assert(j == Jrot.getNumRows());
					
					Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
					float alpha5 = 0.006;
					Matrix JHash = JConcatTranspose * ((ConcatenatedJacobMatrix * JConcatTranspose).inverse());
					Matrix JHashJ = JHash * ConcatenatedJacobMatrix;
					Matrix Intermediate1 = I - JHashJ;
					Vector Intermediate2 = ZeroVector - theta;
					Vector RightTerm = Intermediate1 * Intermediate2;
					Vector LeftTerm = (JHash * ConcatenatedXRDelta) * alpha5;
					Vector deltaTheta = LeftTerm + RightTerm;
					thetahat = deltaTheta + theta;
				}
			}
			else if (stateMachine == MoveTowardsTarget)
			{
				// ----------------------------------------------------------------------------------------
				// Step 5 - Take object, manouvre combined position and orientation control to move towards target GREEN marker.
				//			Orientation control such that BLUE object is oriented w.r.t GREEN marker.
				// ----------------------------------------------------------------------------------------
				Vector rDelta = quatdiff(rhatObject, rhat);
				Vector xdelta = xhat - xhatObject;

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

				assert(j == Jrot.getNumRows());
				
				Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
				float alpha5 = 0.001;
				Matrix JHash = JConcatTranspose * ((ConcatenatedJacobMatrix * JConcatTranspose).inverse());
				Matrix JHashJ = JHash * ConcatenatedJacobMatrix;
				Matrix Intermediate1 = I - JHashJ;
				Vector Intermediate2 = ZeroVector - theta;
				Vector RightTerm = Intermediate1 * Intermediate2;
				Vector LeftTerm = (JHash * ConcatenatedXRDelta) * alpha5;
				Vector deltaTheta = LeftTerm + RightTerm;
				thetahat = deltaTheta + theta;
			}
			else if (stateMachine == OpenGrip)
			{
				// ----------------------------------------------------------------------------------------
				// Step 3 - When close enough to the BLUE object, open claws/set grip and continue combined position and orientation control. Pseudo Inverse Method - Explicit Optimization Criterion, for orientation control
				//			Orientation control remains perpendicular w.r.t BLUE object orientation.
				// ----------------------------------------------------------------------------------------

				Vector rDelta = quatdiff(r, rhatObject);
				Vector xdelta = xhatObject - x;

				// This threshold/constant is based on the size of the blue object
				// and the size of the claws/fingers
				// 
				if (xdelta[0] < 0.036)
				{
					stateMachine = CloseGrip;
				}
				else
				{
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

					assert(j == Jrot.getNumRows());
					
					Matrix JConcatTranspose = ConcatenatedJacobMatrix.transpose();
					float alpha5 = 0.00004;
					Matrix JHash = JConcatTranspose * ((ConcatenatedJacobMatrix * JConcatTranspose).inverse());
					Matrix JHashJ = JHash * ConcatenatedJacobMatrix;
					Matrix Intermediate1 = I - JHashJ;
					Vector Intermediate2 = ZeroVector - theta;
					Vector RightTerm = Intermediate1 * Intermediate2;
					Vector LeftTerm = (JHash * ConcatenatedXRDelta) * alpha5;
					Vector deltaTheta = LeftTerm + RightTerm;
					thetahat = deltaTheta + theta;
				}
			}

			if (stateMachine == OpenGrip)
			{
				// Open the claws of the gripper to be able to encapsulate the object/capsule.
				// 
				setGrip(1, thetahat);
			}
			else if (stateMachine == CloseGrip)
			{		
				// ----------------------------------------------------------------------------------------
				// Step 4 - Close grip, wait till closed to enough extent, based on the target close setting of -0.4
				// ----------------------------------------------------------------------------------------
				// Close the grip to the extent of -0.4 as gripAmount.This number was arrived at
				// as a function of the sizes of the capsule and the claw dimensions.
				// 
				lTarget = 0.4;
				rTarget = -0.4;
				setGrip(-0.4, thetahat);
				stateMachine = ClosingGrip;
			}
			else if (stateMachine == ClosingGrip)
			{
				// If the grip is almost close to the target values (lTarget, rTarget) at this point,
				// Proceed with the next step - Moving towards the target.
				// Else, proceed with ClosingGrip state.
				// 
				if (abs(theta[7] - lTarget) < 0.01 && abs(rTarget - theta[8]) < 0.01)
				{					
					stateMachine = MoveTowardsTarget;
				}
			}

			// set target DOFs to thetahat and advance simulation.
			// 
			mjSetControl(dimtheta, thetahat);
			mjStep2();
		}

		mjDisconnect();
	}
	mjClear();
}