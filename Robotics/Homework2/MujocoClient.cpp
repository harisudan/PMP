// mujoco communication interface
#include "comm.h"

// algebra classes
#include "GMatrix.h"
#include "GVector.h"
#include "Quaternion.h"
#include "windows.h"
#include "conio.h"

typedef GVector<double> Vector;
typedef GMatrix<double> Matrix;
typedef Quaternion<double> Quat;

// useful for neighbor storage and sorting
#include <vector>
#include <algorithm>
#include <forward_list>
#include <queue>
#include <iostream>
#include <fstream>

using namespace std;

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
bool isValidState(const Vector &q, const mjSize &size, bool fSleep)
{
	Vector qvel(size.nv);
	Vector act(size.na);
	qvel.setConstant(0.0);
	act.setConstant(0.0);

	mjSetState(size.nq, size.nv, size.na, 0.0, (Vector&)q, qvel, act);
	mjContact contact = mjGetContacts();

	if (fSleep)
	{
		Sleep(1000);
	}

	return (contact.nc == 0);
}

// test for collision-free path between two body configuration vectors
// return true if there exists a collision-free path between body configuration qa and qb
// return false otherwise
bool isValidPath(const Vector &qa, const Vector &qb, const mjSize &size, bool fSleep, bool fForFinalSimulation)
{
	// number of interpolated steps between qa and qb to test for validity
	static const int numstep = 10;
	static const int dof = 7;

	int i = 0;

	bool isValid = isValidState(qb, size, fSleep);
	assert(isValid);	

	if (!fForFinalSimulation)
	{
		isValid = isValidState(qa, size, fSleep);
		assert(isValid);
	}

	Vector qpos(dof);

	assert(dof == qb.getSize());
	assert(dof == qa.getSize());

	double qaarr[dof] = {0.0};
	double qbarr[dof] = {0.0};

	for (int h = 0; h < dof; h++)
	{
		qaarr[h] = qa[h];
		qbarr[h] = qb[h];
	}

	double interpolateddofvalues[dof][numstep] = { 0 };
	int dofiter = 0;
	double base = 0.0;

	while (dofiter < dof)
	{
		double dofinterpolationstep = (qa[dofiter] - qb[dofiter]) / (numstep + 1);

		// populate ten different in-between points for this particular DOF coordinate.
		// 
		for (int j = 0; j < numstep; j++)
		{
			interpolateddofvalues[dofiter][j] = qb[dofiter] + ((j + 1 ) *dofinterpolationstep);
		}

		dofiter++;
	}

	dofiter = 0;

	while (i < numstep)
	{
		for (int k = 0; k < dof; k++)
		{
			qpos[k] = interpolateddofvalues[k][i];
		}

		if (!isValidState(qpos, size, fSleep))
		{
			return false;
		}

		i++;
	}

	if (fForFinalSimulation)
	{
		isValid = isValidState(qa, size, fSleep);
		assert(isValid);
	}

	return true;
}

// Animate path, by just calling IsValidPath (which already interpolates 10 in-between states and animates them through the call to IsValidState) 
// 
void AnimatePath(const Vector &qa, const Vector &qb, const mjSize &size)
{
	bool isAnimationSuccessful = isValidPath(qa, qb, size, true /*fSleep*/, true /*fForFinalSimulation*/);
	assert(isAnimationSuccessful == true);
}

// forward declarations
// 
class KSample;
typedef KSample* pKSample;
bool GetMinElementFromDijstraQueue(const forward_list<pKSample>& validStatesList, int iTotalSamples, pKSample& pkSampleOut, int sourceIndex);
void ReAdjustEstimateDijkstra(const pKSample& pSample1, pKSample pSample2, int sourceIndex);

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

	pKSample GetPKSample()
	{
		return m_pkSample;
	}

	bool operator<(const QueueEntry& rhs) const
	{
		return m_orderKey < rhs.GetKey();
	}
};

typedef QueueEntry* PQueueEntry;

// A comparator function to compare two entries in the priority queue 
// used for holding K nearest neighbours
// 
struct MyComparator1
{
	bool operator() (PQueueEntry sample1, PQueueEntry sample2)
	{
		// call the overloaded operator
		// 
		return *sample1 < *sample2;
	}
};

// Encapsulates a state as a node in the graph
// 
class KSample
{
	// The number of source states we are dealing with in the problem
	// 
	static const int noOfSource = 3;

	// The state actually represented by this KSample
	// 
	Vector m_sample;

	// The priority queue to carry the K Nearest neighbours for this sample
	// 
	priority_queue<PQueueEntry, vector<PQueueEntry>, MyComparator1> m_KNN;

	// The set of valid transitions (a subset of entries of K Nearest neighbors)
	// populated based on precense of valid 10 interpolated states.
	// 
	forward_list<PQueueEntry> m_validTransitions;

	// The number of valid transitions we found for this sample.
	// 
	int m_iNumberOfValidTransitions;

	// The array carrying the shortest path estimate of this state from each of the three source states.
	// 
	double m_dShortestPathEstimate[noOfSource];

	// The array carrying the backpointer to the previous state used to arrive to this state
	// yielding the above shortest estimate. Array to maintain this information for each of the three source states.
	// 
	pKSample m_pprevSampleWithShortestPathEstimate[noOfSource];

	// A marker to indicate that this sample was previously disjkstra processed (for each of the three sources)
	// 
	bool m_dijktraQueueMarked[noOfSource];

	public:
	KSample(Vector sample)
	{
		m_sample = sample;
		m_iNumberOfValidTransitions = 0;
		m_dShortestPathEstimate[0] = INFINITY;
		m_dShortestPathEstimate[1]= INFINITY;
		m_dShortestPathEstimate[2] = INFINITY;
		m_pprevSampleWithShortestPathEstimate[0] = nullptr;
		m_pprevSampleWithShortestPathEstimate[1] = nullptr;
		m_pprevSampleWithShortestPathEstimate[2] = nullptr;
		m_dijktraQueueMarked[0] = false;
		m_dijktraQueueMarked[1] = false;
		m_dijktraQueueMarked[2] = false;
	}

	Vector GetSampleVector() const
	{
		return m_sample;
	}

	// Function called if this sample is going to be the source state.
	// 
	void SetAsSourceOfShortestPath(int sourceIndex)
	{
		assert(sourceIndex <= noOfSource - 1);
		assert(m_dShortestPathEstimate[sourceIndex] == INFINITY);
		assert(m_dijktraQueueMarked[sourceIndex] == false);
		m_dShortestPathEstimate[sourceIndex] = 0.0;
	}

	double GetShortestPathEstimate(int sourceIndex) const
	{
		assert(sourceIndex <= noOfSource - 1);
		return m_dShortestPathEstimate[sourceIndex];
	}

	// Return the backpointer for a particular source state (based on sourceindex)
	// Used while backtracking and animating the shortest path.
	// 
	void GetPrevSampleInPath(pKSample& pkSample, int sourceIndex)
	{
		assert(sourceIndex <= noOfSource - 1);
		pkSample = m_pprevSampleWithShortestPathEstimate[sourceIndex];
	}

	int GetNumberOfValidTransitions() const
	{
		return m_iNumberOfValidTransitions;
	}

	forward_list<PQueueEntry>& GetValidTransitionsList()
	{
		return m_validTransitions;
	}

	bool IsDijkstraQueueMarked(int sourceIndex) const
	{
		assert(sourceIndex <= noOfSource - 1);
		return m_dijktraQueueMarked[sourceIndex];
	}

	void SetDijkstraQueueMarked(int sourceIndex)
	{
		assert(sourceIndex <= noOfSource - 1);
		assert(!m_dijktraQueueMarked[sourceIndex]);
		m_dijktraQueueMarked[sourceIndex] = true;
	}

	// Update the backpointer and the new shortestpath estimate achieved through that.
	// 
	void SetNewShortestPath(const pKSample& pPrevSampleOnShortestPath, double shortestPathEstimate, int sourceIndex)
	{
		assert(sourceIndex <= noOfSource - 1);
		assert(nullptr != pPrevSampleOnShortestPath);
		m_pprevSampleWithShortestPathEstimate[sourceIndex] = pPrevSampleOnShortestPath;
		m_dShortestPathEstimate[sourceIndex] = shortestPathEstimate;
	}

	// Add a newly generated sample to priority queue, with a max cap on K elements in the queue
	// Retain K closest samples in the queue.
	// 
	void TryAddSampleToPriorityQueue(pKSample pSample, double diffKey, int MaxSizeK)
	{
		assert(nullptr != pSample);
		assert(diffKey > 0.0);

		if (m_KNN.size() >= MaxSizeK)
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

	// Go through each of the K nearest neighbors and find out which of those 
	// can be reached with valid transitions
	// 
	void CheckValidTransitions(const mjSize& size, int MaxSizeK)
	{
		assert(m_KNN.size() == MaxSizeK);

		while (m_KNN.size() > 0)
		{
			PQueueEntry pNeighBourQueueEntry = m_KNN.top();
			assert(nullptr != pNeighBourQueueEntry);

			// Check if there is a valid path to this neighbour
			// 
			if (isValidPath(m_sample, pNeighBourQueueEntry->GetPKSample()->GetSampleVector(), size, false /*fSleep*/, false /*fForFinalSimulation*/))
			{
				m_iNumberOfValidTransitions++;
				m_validTransitions.emplace_front(new QueueEntry(*pNeighBourQueueEntry));
			}
			m_KNN.pop();
		}

		assert(m_iNumberOfValidTransitions <= MaxSizeK);
	}

	// Function to unit test that all transitions in m_validTransitions are indeed valid.
	// 
	void TestValidTransitions(const mjSize& size)
	{
		auto it = m_validTransitions.begin();

		for (int i = 0; i < m_iNumberOfValidTransitions; i++)
		{
			PQueueEntry pQueueEntryValidTransitionTarget = (PQueueEntry)(*it++);
			assert(nullptr != pQueueEntryValidTransitionTarget);

			assert(isValidPath(m_sample, pQueueEntryValidTransitionTarget->GetPKSample()->GetSampleVector(), size, false /*fSleep*/, false /*fForFinalSimulation*/));
		}
	}
};

// Execute Dijskstra's algorithm to find the shortest paths between pkSource and any vertex.
// This would also compute the shortest path to pKGoal.
// 
void FindShortestPath(int sourceIndex, pKSample pkSource, forward_list<pKSample>& validStates, int iTotalSamples)
{
	// Test the valid transitions
	// 
	pKSample pk = nullptr;
	int y = 0, j = 0;
	auto it = validStates.begin();

	// Asserting the condition of the shortestPathEstimate variable and 
	// marking the source's shortestPathEstimate as ZERO.
	// 
	for (j = 0; j < iTotalSamples; j++)
	{
		pk = (pKSample)(*it++);
		assert(nullptr != pk);
		
		assert(pk->GetShortestPathEstimate(sourceIndex) == INFINITY);

		// Mark the shortestpath estimate for the source state as zero.
		// 
		if (pk == pkSource)
		{
			y++;
			pk->SetAsSourceOfShortestPath(sourceIndex);
		}
	}

	// Only one source - just a sanity to guard against any bad pointer management.
	// 
	assert(y == 1);

	pKSample pkSampleDijkstraIterator = nullptr;

	// No:of states in the graph that are reachable, just to keep track.
	// 
	int noOfReachableStates = 0;

	// While we have connected nodes in the graph which have not been smoothed before,
	// retrieve each vertex based on whosoever has the min shortest path estimate (simulating a min priority queue)
	// and ReAdjust (update shortest path if reaching the neighbour is shorter through this node) all the links that the vertex has.
	// 
	while (GetMinElementFromDijstraQueue(validStates, iTotalSamples, pkSampleDijkstraIterator, sourceIndex))
	{
		noOfReachableStates++;
		assert(nullptr != pkSampleDijkstraIterator);

		forward_list<PQueueEntry> validTransitionsForSample = pkSampleDijkstraIterator->GetValidTransitionsList();

		auto it = validTransitionsForSample.begin();
		for (int p = 0; p < pkSampleDijkstraIterator->GetNumberOfValidTransitions() && it != validTransitionsForSample.end(); p++)
		{
			PQueueEntry pkNeighbour = (PQueueEntry)(*it);
			assert(nullptr != pkNeighbour);

			// Try to see if reaching the pkNeighbour through pkSampleDijkstraIterator will yield in a shorter path
			// If so, update the shortest path estimate and the back pointer to point to this parent (for backtracking)
			// 
			ReAdjustEstimateDijkstra(pkSampleDijkstraIterator, pkNeighbour->GetPKSample(), sourceIndex);
			it++;
		}

		// Mark that we have processed this vertex once and so will not process it anymore.
		// 
		pkSampleDijkstraIterator->SetDijkstraQueueMarked(sourceIndex);
	}

	// Unit Test the disjkstra.
	// Assert some of the expected facts to make sure all objects and data structures are in expected states.
	// 
	pKSample prevSampleInPath = nullptr;
	it = validStates.begin();
	int noOfUnreachableStates = 0;
	for (j = 0; j < iTotalSamples; j++)
	{
		pk = (pKSample)(*it++);
		assert(nullptr != pk);

		pk->GetPrevSampleInPath(prevSampleInPath, sourceIndex);

		// All vertices which are still marked with shortestpathestimate as INFINITY are unreachable nodes in the graph.
		// 
		if (pk->GetShortestPathEstimate(sourceIndex) != INFINITY)
		{
			assert(pk->IsDijkstraQueueMarked(sourceIndex) == true);
			assert(prevSampleInPath != nullptr || (pk->GetShortestPathEstimate(sourceIndex) == 0.0 && pk == pkSource));
		}
		else
		{
			// these are unreachable states from the other validStates.
			// 
			assert(pk->IsDijkstraQueueMarked(sourceIndex) == false);
			assert(prevSampleInPath == nullptr);
			noOfUnreachableStates++;
		}
	}
}

// Given the source state and the goal state, animate the path going from the source to the target
// Use the backpointers stored in each pKSample to get to the pkSource. In-between each sample on the way,
// use the IsValidPath function to interpolate 10 valid states in-between those two states and animate them.
// 
double AnimateShortestPaths(int sourceIndex, pKSample pkSource, pKSample pKGoal, const mjSize& size)
{
	// get input from user first so they are ready for seeing the animated reconstruction of shortest paths
	//	
	cout << "Press enter key to do shortest path animation for source: " << sourceIndex << "\n";
	char c = getchar();

	// Test the valid transitions
	// 
	pKSample pk = nullptr;
	pKSample pkCurrentBackIterator = nullptr;
	int y = 0, j = 0;
	double costOfPath = 0.0;

	assert(nullptr != pKGoal);
	assert(nullptr != pkSource);
	assert(sourceIndex <= 2);

	pkCurrentBackIterator = pKGoal;
	forward_list<pKSample> forwardPathOrder;
	int noOfStatesInPath = 0;

	// Go through the back pointers from pKGoal and add them in to forwardPathOrder list.
	// So the list is populated such a way that the state ending up at the front of the list is the source.
	// and the state at the end of the list is the goal.
	// 
	while (pkCurrentBackIterator != pkSource)
	{
		forwardPathOrder.emplace_front(pkCurrentBackIterator);
		noOfStatesInPath++;
		pkCurrentBackIterator->GetPrevSampleInPath(pk, sourceIndex);
		assert(nullptr != pk || (pk == pkSource && (pk->GetShortestPathEstimate(sourceIndex) == 0.0)));
		pkCurrentBackIterator = pk;
	}

	// Place the source at the starting of the list.
	// 
	assert(pkCurrentBackIterator == pkSource);
	forwardPathOrder.emplace_front(pkSource);
	noOfStatesInPath++;

	// Now iterate forward starting from the front (the source state) and animate each state (with 10 interpolated values between each state).
	// 
	auto it = forwardPathOrder.begin();

	pKSample pkPrevious = nullptr;

	for (j = 0; j < noOfStatesInPath; j++, it++)
	{
		pk = (pKSample)(*it);
		assert(nullptr != pk);

		if (nullptr != pkPrevious)
		{
			// this interpolates 10 points inbetween and simulates it
			// also asserts that isValidPath returns true
			// 
			AnimatePath(pk->GetSampleVector(), pkPrevious->GetSampleVector(), size);
		}
		else
		{
			// source state (no previous state)
			// 
			assert(j == 0 && pk == pkSource);
			Vector qvel(size.nv);
			Vector act(size.na);
			qvel.setConstant(0.0);
			act.setConstant(0.0);

			mjSetState(size.nq, size.nv, size.na, 0.0, (Vector&)(pk->GetSampleVector()), qvel, act);
			mjContact contact = mjGetContacts();
			Sleep(500);

			assert(contact.nc == 0);
		}

		pkPrevious = pk;
	}

	assert(pk == pKGoal);
	assert(pkPrevious == pKGoal);

	// Return the shortest path estimate for this source (represented by the sourceindex) from the goal.
	// 
	return pKGoal->GetShortestPathEstimate(sourceIndex);
}

// Helper function used in Djkstra's update step.
// Try to see if reaching pSample2 from pSample1 is closer than the shortest path estimate we already have for pSample2.
// If yes, update the new shortest path estimate for reaching through pSample1 and also 
// update pSample2's back pointer to pSample1, to assist in correct back tracking.
// 
void ReAdjustEstimateDijkstra(const pKSample& pSample1, pKSample pSample2, int sourceIndex)
{
	assert(nullptr != pSample1);
	assert(nullptr != pSample2);

	Vector pSample1Vector = pSample1->GetSampleVector();
	Vector pSample2Vector = pSample2->GetSampleVector();

	// Take the diff between the two vectors and take the norm of the vector.
	// 
	Vector diff = pSample1Vector - pSample2Vector;

	double distanceBetween1and2 = diff.length();
	assert(distanceBetween1and2 > 0.0);

	// If the estimate achieved throughpSample1 is lower than the existing shortest path estimate for pSample2, update the path and estimate.
	// 
	if (pSample1->GetShortestPathEstimate(sourceIndex) + distanceBetween1and2 < pSample2->GetShortestPathEstimate(sourceIndex))
	{
		pSample2->SetNewShortestPath(pSample1, pSample1->GetShortestPathEstimate(sourceIndex) + distanceBetween1and2, sourceIndex);
	}
}

// This is the function simulating a min priority queue based on the shortestpath estimates of every state (node in a graph)
// Simply go through the existing list in validStatesList, whosoever has the shortest estimate in the list and 
// whosoever was not dijkstra processed before, return that state.
// 
bool GetMinElementFromDijstraQueue(const forward_list<pKSample>& validStatesList, int iTotalSamples, pKSample& pkSampleOut, int sourceIndex)
{	
	auto it = validStatesList.begin();
	double minShortestPathEstimateFound = INFINITY;

	for (int i = 0; i < iTotalSamples && it != validStatesList.end(); i++, it++)
	{
		pKSample pk = pKSample(*it);
		assert(nullptr != pk);

		if (pk->GetShortestPathEstimate(sourceIndex) < minShortestPathEstimateFound && !pk->IsDijkstraQueueMarked(sourceIndex))
		{
			pkSampleOut = (*it);
			minShortestPathEstimateFound = pk->GetShortestPathEstimate(sourceIndex);
		}
	}

	// don't process un-reachable nodes
	// 
	return minShortestPathEstimateFound != INFINITY;
}

void main(int argc, char* argv[])
{
	// three initial poses vectors (must be included in the list of samples)
	static const Vector qinit[3] = { Vector(NUM_JOINT, keyinit[0]), 
		Vector(NUM_JOINT, keyinit[1]), Vector(NUM_JOINT, keyinit[2]) };

	// goal pose vector (must be included in the list of samples)
	static const Vector qgoal = Vector(NUM_JOINT, keygoal);

	static const int N = argc > 1 ? atoi(argv[1]) : 1000;
	static const int K = argc > 2 ? atoi(argv[2]) : 30;
	
	// + 4 to account for 3 init states and 1 goal state.
	// 
	static const int iTotalSamples = N + 4;

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
		// 
		forward_list<pKSample> validStates(iTotalSamples + 4);
		Vector qpos;
		pKSample pkSample = nullptr;
		int i = 0, j = 0;

		// Generate N samples and add them to validStates if they are valid.
		// 
		while (i < N)
		{
			qpos = generateSample();
			
			// Check if its collision free.
			// 
			if (isValidState(qpos, size, false))
			{
				pkSample = new KSample(qpos);
				validStates.emplace_front(pkSample);
				i++;
			}
		}
		
		assert(i == N);
		pKSample pkqinit0 = new KSample(qinit[0]);
		pKSample pkqinit1 = new KSample(qinit[1]);
		pKSample pkqinit2 = new KSample(qinit[2]);
		pKSample pkqgoal = new KSample(qgoal);

		// Add the three init and 1 goal state to the list of validStates.
		// 
		validStates.emplace_front(pkqinit0);
		validStates.emplace_front(pkqinit1);
		validStates.emplace_front(pkqinit2);
		validStates.emplace_front(pkqgoal);

		// point the iterator the beginning of the forwardlist.
		// 
		auto it = validStates.begin();
		pKSample pk = nullptr;
		
		for (j = 0; j < iTotalSamples; j++)
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

		// loop through all samples to generate K closest neighbours for each sample
		// 
		for (j = 0; j < iTotalSamples; j++)
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

				// Try to insert the element with the diffLength distance in to priority queue.
				// The function will take care to retain only K nearest neighbors
				// 
				pkSampleCurrent->TryAddSampleToPriorityQueue(pk, diffLength, K);

				// advance the it pointer
				//
				it++;
			} while (++k < iTotalSamples);

			// Calculate the nearest neighbours of the next sample.
			// 
			currentSample++;
		}

		// find the valid transitions for each sample from the set of K nearest neighbours
		//
		it = validStates.begin();
		for (j = 0; j < iTotalSamples; j++)
		{
			pk = (pKSample)(*it++);
			assert(nullptr != pk);

			// this will loop through each of K nearest neighbours and populate the validTransitions list
			// 
			pk->CheckValidTransitions(size, K);
		}

		// Test the valid transitions
		// 
		it = validStates.begin();
		for (j = 0; j < iTotalSamples; j++)
		{
			pk = (pKSample)(*it++);
			assert(nullptr != pk);

			pk->TestValidTransitions(size);
		}
		
		// find three shortest paths among the available samples from 
		// each of three qinit vectors to qgoal vector (using nearest
		// neighbors above)
		// 
		FindShortestPath(0, pkqinit0, validStates, iTotalSamples);
		FindShortestPath(1, pkqinit1, validStates, iTotalSamples);
		FindShortestPath(2, pkqinit2, validStates, iTotalSamples);
		

		// animate three shortest paths found above by interpolating their
		// positions and setting this interpolated state by mjSetState()
		// 
		double costFromSource[3] = { 0.0 };
		costFromSource[0] = AnimateShortestPaths(0, pkqinit0, pkqgoal, size);
		costFromSource[1] = AnimateShortestPaths(1, pkqinit1, pkqgoal, size);
		costFromSource[2] = AnimateShortestPaths(2, pkqinit2, pkqgoal, size);

		// done!
		// write down the results (shortest path cost for each of three sources to the goal)
		// 
		if (argc > 3)
		{
			ofstream myfile;
			myfile.open(argv[3], ios_base::app);
			myfile << "SampleSize: " << N << " K: "<< K << "\n";
			myfile << "Shortest Path Estimate For Each Source, for the above parameters.\n";
			myfile << "Source: Qinit0; ShortestPathDistanceEstimate: " << costFromSource[0] << "\n";
			myfile << "Source: Qinit1; ShortestPathDistanceEstimate: " << costFromSource[1] << "\n";
			myfile << "Source: Qinit2; ShortestPathDistanceEstimate: " << costFromSource[2] << "\n";
			myfile << "\n\n";
			myfile.close();
		}
		mjDisconnect();
	}
	mjClear();
}
