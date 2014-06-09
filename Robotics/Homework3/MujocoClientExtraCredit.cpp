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

// paddle geom ID
static const int PADDLE_GEOM = 4;
// number of trial keyframes in the model
static const int NUM_TRIAL = 10;
// simulation timestep size
static const double dt = 1e-3;

// model dynamics function
// given x{k-1} returns x{k} = f(x{k-1}) and its Jacobian J_f(x{k-1})
// J_f is |x|*|x| matrix
void f(const Vector &xprev, Vector *x, Matrix *Jfxprev)
{
	// acceleration due to gravity
	static const double a[3] = { 0.0, 0.0, -2.0 };

	// generate linear dynamics matrices
	Matrix A(6, 6);
	for (int r = 0; r<3; r++)
	for (int c = 0; c<3; c++)
	{
		A(0 + r, 0 + c) = (r == c) ? 1.0 : 0.0;
		A(3 + r, 3 + c) = (r == c) ? 1.0 : 0.0;
		A(0 + r, 3 + c) = (r == c) ? dt : 0.0;
		A(3 + r, 0 + c) = 0.0;
	}
	Vector b(6);
	for (int r = 0; r<3; r++)
	{
		b[0 + r] = 0.0;
		b[3 + r] = dt * a[r];
	}

	if (x)
	{
		*x = A * xprev + b;
	}
	if (Jfxprev)
	{
		*Jfxprev = A;
	}
}

// observation function
// given x{k} return z{k} = h(x{k}) and its Jacobian J_h(x{k})
// J_h is |z|*|x| matrix 
void h(const Vector &x, Vector *z, Matrix *Jhx)
{
	double k1 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	double k2 = x[0] * x[0] + x[1] * x[1];

	if (z)
	{
		*z = Vector(3);
		(*z)[0] = atan2(x[1], x[0]);
		(*z)[1] = atan2(x[2], sqrt(k2));
		(*z)[2] = sqrt(k1);
	}
	if (Jhx)
	{
		*Jhx = Matrix(3, 6);
		Jhx->setConstant(0.0);

		// dz[0]/dx
		(*Jhx)(0, 0) = -x[1] / k2;
		(*Jhx)(0, 1) = +x[0] / k2;
		(*Jhx)(0, 2) = 0.0;

		// dz[1]/dx
		(*Jhx)(1, 0) = -x[0] * x[2] / (sqrt(k2) * k1);
		(*Jhx)(1, 1) = -x[1] * x[2] / (sqrt(k2) * k1);
		(*Jhx)(1, 2) = sqrt(k2) / k1;

		// dz[2]/dx
		(*Jhx)(2, 0) = x[0] / sqrt(k1);
		(*Jhx)(2, 1) = x[1] / sqrt(k1);
		(*Jhx)(2, 2) = x[2] / sqrt(k1);
	}
}

void main(void)
{
	// dynamics noise covariance matrix
	// dynamics are deterministic in this homework, so Q = 0
	Matrix Q(6, 6);
	Q.setConstant(0.0);

	// observation noise covariance matrix
	// change the diagonal entries of this if you change noise field in xml model
	Matrix R(3, 3);
	R.setIdentity();
	R(0, 0) = 0.05;
	R(1, 1) = 0.05;
	R(2, 2) = 0.30;

	// connect to mujoco server
	mjInit();
	mjConnect(10000);

	// load hand model
	if (mjIsConnected())
	{
		mjLoad("estimator.xml");
		Sleep(1000);  // wait until the load is complete
	}

	if (mjIsConnected() && mjIsModel())
	{
		mjSetMode(2);

		// loops over a series of simulation trials (state is reset between trials)
		int numHit = 0;
		int numMiss = 0;

		// identity matrix for re-use
		//
		Matrix iden(6, 6);
		iden.setIdentity();

		for (int trial = 0; trial<NUM_TRIAL; trial++)
		{
			// reset to trial's keyframe
			mjReset(trial);

			// estimate of initial state
			// set xa to something non-informative, but non-zero to prevent singularities
			Vector xa(6);
			xa.setConstant(0.0); xa[1] = 1.0;

			// noise covarience in initial state estimate
			Matrix P(6, 6);

			// TODO: play with this parameter to see how it affects results
			double Pdiag = 1e-0;
			P.setIdentity();
			P = P * Pdiag;

			Matrix Pk = P;

			// run simulation until ball hits something
			for (int k = 0;; k++)
			{
				// noisy sensor measurement
				Vector z(3, mjGetSensor());

				double blah = z[0];
				double blah1 = z[1];
				double blah2 = z[2];

				// TODO:
				// perform your EFK predictor/corrector updates here
				// you should have xa and P updated after this

				// Model forecast step
				// 
				Vector xkf(6);
				Matrix jfxk(6, 6);
				f(xa, &xkf, &jfxk);

				Matrix Pkf(6, 6);
				Pkf.setConstant(0.0);
				Pkf = (jfxk * Pk * jfxk.transpose()) + Q;

				// Data assimilation step
				// 
				Matrix kalmanGain_i(6, 6);
				Vector hxk_i(3);
				Matrix jhxk_i(6, 6);
				Vector xk_i = xkf;
				Vector xk_i_previous = xk_i;
				double diff = 0.0;

				for (int i = 0; i < 40; i++)
				{
					h(xk_i, &hxk_i, &jhxk_i);

					// calculate kalman gain
					// 
					kalmanGain_i = Pkf * jhxk_i.transpose() * ((jhxk_i * Pkf * jhxk_i.transpose() + R).inverse());

					// Update xk_i
					// 
					xk_i = xkf + kalmanGain_i * (z - hxk_i);

					// Calculate diff (change) from previos xk_i_previous
					// Stop the loop if the diff (change) is less
					//
					diff = (xk_i - xk_i_previous).length();

					// A threshold of 1e-4 used to stop on further iterations
					// 
					if (diff < 0.0001)
					{
						break;
					}

					// update xk_i_previous
					//
					xk_i_previous = xk_i;
				}

				Pk = (iden - (kalmanGain_i * jhxk_i)) * Pkf;
				xa = xk_i;

				// set estimator location for visualization.
				// 
				mjSetEstimator(xa);

				Vector xtarget(2);

				// TODO:
				// fill xtarget vector here with desired horizonal (entry 0) and
				// vectical (entry 1) position of the paddle
				// set controls.

				// Calculate timeToImpact by dividing Py / vy (there is no acceleration across y)
				//
				double timeToImpact = xa[1] / xa[4];

				if (timeToImpact < 0)
				{
					timeToImpact = -timeToImpact;
				}

				timeToImpact = timeToImpact / 1000;

				double blah3 = xa[0];
				double blah4 = xa[1];
				double blah5 = xa[2];
				double blah6 = xa[3];
				double blah7 = xa[4];
				double blah8 = xa[5];

				// Calculate the height at impact (the z coordinate where the ball is estimated to hit the wall)
				// by using the timeToImpact, Pz, vz and accelaration due to gravity (-2 in z direction)
				//
				double zHeightAtImpact = xa[2] + (xa[5] * timeToImpact) + (0.5 * -2 * (timeToImpact * timeToImpact));

				// Do similar calculation to obtain x coordinate at impact - but the acceleration across x is zero
				// 
				double xCoordinateAtImpact = xa[0] + (xa[3] * timeToImpact);

				// Now set xtarget position of paddle to the estimated impact coordinates of xCoordinateAtImpact, zHeightAtImpact
				// 
				xtarget[0] = xCoordinateAtImpact;
				xtarget[1] = zHeightAtImpact;

				mjSetControl(2, xtarget);

				// step simulation forward
				mjStep();
				Sleep(1);

				// check for paddle of wall contact to determine scoring
				// do not modify this
				mjContact contact = mjGetContacts();
				if (contact.nc > 0)
				{
					printf("trial %02d: %02d ", trial, int(contact.geom1[0]));
					if (contact.geom1[0] == PADDLE_GEOM)
					{
						printf("hit\n");
						numHit++;
					}
					else
					{
						printf("miss\n");
						numMiss++;
					}

					// run the simulation for one more second to show post-contact result
					for (int kpost = 0; kpost<int(1.0 / dt); kpost++)
					{
						mjStep();
						Sleep(1);
					}

					break;
				}
			}
		}
		// we will use this to tally your results. Do not remove.
		printf("hit: %02d miss: %02d\n", numHit, numMiss);

		mjDisconnect();
	}
	mjClear();
}