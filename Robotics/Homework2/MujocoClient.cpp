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

// useful for neighbor storage and sorting
#include <vector>
#include <algorithm>
#include <forward_list>
#include <queue>

#ifndef INFINITY
// you can use this constant to denote infinitely large values
static const double INFINITY = DBL_MAX;
#endif
// hardcoded number of model joints
static const int NUM_JOINT = 7;
// joint limits for the model
static const double jointmin[NUM_JOINT] = {-40, -50, -50, -60, -90, -30, -60};
static const double jointmax[NUM_JOINT] = {+40, +50, +50, +60, +90, +30, +60};
// initial poses
static const double keyinit[3][NUM_JOINT] = 
{
	{0.35, -0.22, 0.83, -0.42, 0.84, 0, 0.35},
	{0.7, 0.87, 0.87, 1, -0.38, 0.17, 0},
	{-0.64, -0.21, -0.47, -0.1, -0.12, 0, 0}
};
// goal pose
static const double keygoal[NUM_JOINT] = {0, 0, 0, 0, 0, 0, 0};


// generate a position vector for a random arm configuration (within joint limits)
// returns a qpos vector suitable for use in mjSetState()
Vector generateSample()
{
	static const double DEG_TO_RAD = 0.017453292519943295769236907684886;

	Vector qpos(NUM_JOINT);
	for(int i=0; i<NUM_JOINT; i++)
	{
		qpos[i] = double(rand()) / double(RAND_MAX);
		qpos[i] = (jointmax[i] - jointmin[i]) * qpos[i] + jointmin[i];
		qpos[i] *= DEG_TO_RAD;
	}
	return qpos;
}

// test for a collision-free body configuration vector
// return true if body configuration q is collision-free
// return false otherwise
// do not change this function
bool isValidState(const Vector &q, const mjSize &size)
{
	Vector qvel(size.nv);
	Vector act(size.na);
	qvel.setConstant(0.0);
	act.setConstant(0.0);

	mjSetState(size.nq, size.nv, size.na, 0.0, (Vector&)q, qvel, act);
	mjContact contact = mjGetContacts();
	return (contact.nc == 0);
}

// test for collision-free path between two body configuration vectors
// return true if there exists a collision-free path between body configuration qa and qb
// return false otherwise
bool isValidPath(const Vector &qa, const Vector &qb, const mjSize &size)
{
	// number of interpolated steps between qa and qb to test for validity
	static const int numstep = 10;

	// your implementation goes here
	return true;
}

class KSample;
typedef KSample* pKSample;

class QueueEntry
{
	pKSample m_pkSample; 
	double m_orderKey;

public:
	QueueEntry(pKSample pkSample, double orderKey)
	{
		m_pkSample = pkSample;
		m_orderKey = orderKey;
	}

	double GetKey() const
	{
		return m_orderKey;
	}

	bool operator<(const QueueEntry& rhs) const
	{
		return m_orderKey < rhs.GetKey();
	}
};

typedef QueueEntry* PQueueEntry;

struct MyComparator
{
	bool operator() (PQueueEntry sample1, PQueueEntry sample2)
	{
		return *sample1 < *sample2;
	}
};

class KSample
{
	static const int K = 20;
	Vector m_sample;
	std::priority_queue<PQueueEntry, std::vector<PQueueEntry>, MyComparator> m_KNN;

	public:
	KSample(Vector sample)
	{
		m_sample = sample;
	}

	Vector GetSampleVector()
	{
		return m_sample;
	}

	void TryAddSampleToPriorityQueue(pKSample pSample, double diffKey)
	{
		if (m_KNN.size() >= K)
		{
			// pop the sample farthest from this sample.
			// 
			if (m_KNN.top()->GetKey() > diffKey)
			{
				m_KNN.pop();
				m_KNN.emplace(new QueueEntry(pSample, diffKey));
			}
		}
		else
		{
			m_KNN.emplace(new QueueEntry(pSample, diffKey));
		}		
	}
};

void main(void)
{
	// three initial poses vectors (must be included in the list of samples)
	static const Vector qinit[3] = { Vector(NUM_JOINT, keyinit[0]), 
		Vector(NUM_JOINT, keyinit[1]), Vector(NUM_JOINT, keyinit[2]) };

	// goal pose vector (must be included in the list of samples)
	static const Vector qgoal = Vector(NUM_JOINT, keygoal);

	static const int N = 1001;
	static const int K = 10;
	static const int nSamples = 100;
	static const int nTotalSamples = nSamples + 4;

	// connect to mujoco server
	mjInit();
	mjConnect(10000);

	// load hand model
	if(mjIsConnected())
	{
		mjLoad("hand-HW2.xml");
		Sleep(1000);  // wait till the load is complete
	}

	if(mjIsConnected() && mjIsModel())
	{
		mjSetMode(1);
		mjReset();

		// size containts model dimensions
		mjSize size = mjGetSize();

		// your implemention:

		// generate a collection of valid samples here
		// make sure qinit and qgoal are included in this collection

		std::forward_list<pKSample> validStates(nTotalSamples + 4);
		Vector qpos;
		pKSample pkSample = nullptr;
		int i = 0;

		while (i < nSamples)
		{
			qpos = generateSample();
			
			// Check if its collision free.
			// 
			if (isValidState(qpos, size))
			{
				pkSample = new KSample(qpos);
				validStates.emplace_front(pkSample);
				i++;
			}
		}
		
		assert(i == nSamples);
		pKSample pkqinit0 = new KSample(qinit[0]);
		pKSample pkqinit1 = new KSample(qinit[1]);
		pKSample pkqinit2 = new KSample(qinit[2]);
		pKSample pkqgoal = new KSample(qgoal);

		// Add the three init and 1 goal state.
		// 
		validStates.emplace_front(pkqinit0);
		validStates.emplace_front(pkqinit1);
		validStates.emplace_front(pkqinit2);
		validStates.emplace_front(pkqgoal);

		// point the iterator the beginning of the forwardlist.
		// 
		auto it = validStates.begin();
		pKSample pk = nullptr;
		
		for (int j = 0; j < nTotalSamples; j++)
		{
			pk = (pKSample)(*it++);
			assert(nullptr != pk);

			qpos = pk->GetSampleVector();
			
			// example: here's how to set this qpos in mujoco server.
			// 
			Vector qvel(size.nv);
			Vector act(size.na);
			qvel.setConstant(0.0);
			act.setConstant(0.0);

			mjSetState(size.nq, size.nv, size.na, 0.0, qpos, qvel, act);
			Sleep(100);
		}
							
		// once valid, collision-free samples are generated, create a set of
		// nearest neighbors for each sample.
		//
		Vector qposCurrent;
		Vector diff;
		double diffLength = 0;
		pKSample pkSampleCurrent = nullptr;
		auto currentSample = validStates.begin();
		int k = 0;

		for (int j = 0; j < nTotalSamples; j++)
		{
			k = 0;
			pkSampleCurrent = (pKSample)(*currentSample);
			assert(nullptr != pkSampleCurrent);

			qposCurrent = pkSampleCurrent->GetSampleVector();

			it = validStates.begin();
			do
			{				
				// Don't compare the element against itself, move ahead.
				// 
				if (it == currentSample)
				{
					it++;
					continue;
				}

				pk = (pKSample)(*it);
				assert(nullptr != pk);

				qpos = pk->GetSampleVector();

				// Calcular difference between the two (qposCurrent and qpos).
				//
				diff = qposCurrent - qpos;
				diffLength = diff.length();

				// Try to insert the element with the diffLength in to priority heap.
				// 
				pkSampleCurrent->TryAddSampleToPriorityQueue(pk, diffLength);

				// advance the it pointer
				//
				it++;
			} while (++k < nTotalSamples);

			// Calculate the nearest neighbours of the next sample.
			// 
			currentSample++;
		}
		
		// find three shortest paths among the available samples from 
		// each of three qinit vectors to qgoal vector (using nearest
		// neighbors above)

		// animate three shortest paths found above by interpolating their
		// positions and setting this interpolated state by mjSetState()

		// done!

		mjDisconnect();
	}
	mjClear();
}
