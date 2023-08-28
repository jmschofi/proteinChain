#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <string>
#include <vector>
#include <iterator>
#include <map>
#include <sys/time.h>
#include <iomanip>
#include <random>


#include <ccgsl/function_scl.hpp>
#include <ccgsl/integration.hpp>

//#include <omp.h>
#include "vecmat3.h"
#include "rndgen.h"
#include "report.h"
#include "inifile.hxx"
#include "binfile.hpp"
#include "myinline.h"
#include "stopwatch.h"
#include "state.h"
#include "fptClass.h"


#include "smoothpdfSpline.h"
#include "correlations.h"
#include "histogram.h"
#include "fpt.h"
#include "event.h"
#include "queue.h"
#include "uniqueBond.h"
#include "internalState.h"
#include "configurationPool.h"


template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ""));
  }
  return out;
}

static bool sortByLevel(State &lhs, State &rhs) {

  if (lhs.getLevel() == rhs.getLevel() ) return lhs.getDecimalIndex() < rhs.getDecimalIndex();

  else return lhs.getLevel() < rhs.getLevel();
}


enum EventType { MIN_COLLISION=0, MAX_COLLISION, BOND_COLLISION, CONFIGMEASUREMENT, STOP, SYNCHRONIZE, OUTPUT, DISTANCEMEASUREMENT, ACFMEASUREMENT};
using namespace std;

//typedef std::vector<Vector> Configuration; // set of positions of all atoms

class ProteinChain {
 private:
  std::default_random_engine generator;
  int numAtoms;
  int numInteractions;
  int numAttractions;
  
  std::vector<Vector> pos;
  std::vector<Vector> posTrial;
  std::vector<Vector> posSave;
  std::vector<Vector> vel;
  std::vector<Vector> velSave;
 
  int maxBonds;
  const static int numberBits = 26;
  long int saveState, previousState, currentState; // number assigned to current configuration
  int zeroState; // state used for zero of entropy
  std::vector<bool> levelTarget; // states for distance data
  std::vector<int> stateConnected; // array to indicate count of transitions from previous state from unbonded target can reach next level
  int saveIter;
  double startTime;


  //double* collTimes;
  std::vector<double> localTime; // holds time of local clock for each atom
  double localTimer; // old variable

  std::vector<int> interactionNumber;  // index counter for hybrid queue
  string interactionFilename;

  int numBonds, numStates;
  int numBondsFixed, numStatesFixed;
  int numAllowedStates;

  std::vector<int> bondIndex1;
  std::vector<int> bondIndex2;
  std::vector<bool> isBondedAtom;
  std::vector<bool> bondExists;
  std::vector<bool> allowedState;
  std::vector<State> states;
  std::vector<ProblemPair> problemStates; // states to re-run with staircase
  std::vector<int> stateIndexToState; // map from index to sorted state index

  std::vector<int> bondIndex1_save;
  std::vector<int> bondIndex2_save;
  std::vector<bool> isBondedAtom_save;
  std::vector<bool> bondExists_save;
  std::vector<bool> allowedState_save;
  std::vector<State> states_save;
  std::vector<int> stateIndexToState_save; // map from index to sorted state index
 
  bool bondsInitialized;
  std::vector<long int> configurationCount;
  std::vector<double> bondDistance;
  std::vector<long int> configurationCount_save;
  std::vector<double> bondDistance_save;


  int numBondEvents; // number of config changes in dynamics
  int numBondMC;

  std::map<int,int> localHist;
  //int* index1;
  //int* index2;

  int collider1,collider2;

  double kineticEnergy;
  double potEnergy, totEnergy, initialEnergy, kineticEnergy1;
  double potEnergy1, totEnergy1 ; // to verify the calculation of the energy

  struct timeval tv;
  struct timezone tz;
  int ret0, ret1;
  double sec0, microsec0;
  double sec1, microsec1;
  double elapsedTime;
  double picoTime;
  int counterPicoTime;
  double simTime;

  int numIterations;
  int numMC;
  double timeMD;
  double basicTimeMD;
  double staircaseTime;
  double isolationTime;

  double outputTime, outputInterval, boxLength;
  double syncInterval;
  double maxBondLengthSq1,maxBondLengthSq2;
  double minBondLengthSq_hc, minBondLengthSq1, minBondLengthSq2;
  double innerBondDistance;

  bool doneScanIteration;

  double beta;
  double mass;

  double timeNow;
  double globalTime;

  EventType collisionType;

  double eventTime, totalTime;

  //    Old Queue information that is not used
  double* eventTimes;
  EventType* eventTypes;
  Event *evTree; // The nodes of the tree
  Event evSchedule;

  int currentType;

  bool initialConfigFlag, movieFlag, movieInitialized, acfFlag, redrawVelocityFlag;
  bool configFlag;
  double vij_rel_av, acfInterval, acfTotalInterval, acfCorrelationTime;
  int acfPoints;
  Correlation<double, double> *acf;

  SmoothPdfSpline** spdf;

  int numAv;
  double Eav, PEav, KEav;

  double checkpointInterval, checkpointTime, redrawVelocityTime, redrawVelocityInterval;
  double runTime;
  double configInterval;
  int configMCInterval;
  int numCycles; // number of iterations to update entropy
  int  configsPerState;

  double timeError, errorInterval; // for testing purposes
  bool errorFlag;
  int errorCount;
  double movieInterval, movieTime;

  bool equilibrateFlag;
  double equilibrationTime;

  //  hybrid scheduler
  int nlists;
  int* linearLists;
  int currentIndex;
  double baseIndex;
  int numInTree, numInLists, insertCounter, scale;
  //
  int poolSize, evIdA, evIdB;
  //
  //  New queue stuff
  //
  HybridQueue* hybridQ;

  Inifile settings;
  int randomSeed;

  bool dataFlag;
  double dataInterval, dataTime;
  int numDataPoints;
  binfile dataFile;
  std::vector<float> distances;
  std::vector<std::vector<float>> distanceMatrix;
  int maxNumDistances;
  
  bool inner;
  //std::vector< std::vector<bool> > innerDone;// record whether inner fpt has been done for pair
  //std::vector< std::vector<bool> > outerDone;
  //std::vector< std::vector<double> > innerFPT;
  //std::vector< std::vector<double> > outerFPT;  // store values to determine optimal pairs and staircases
  
  std::vector< std::vector<FptData> > innerFPT;
  std::vector< std::vector<FptData> > outerFPT;


  bool biasFlag;
  string biasFilename;
  std::vector<double> Sbias;
  std::vector<double> Sbias_save;
  std::vector<double> dS;
  double pcrit;
  double entropyScaler,maxEntropyScaler,minEntropyScaler;
  std::vector<std::vector<double>> entropyMatrix;
  int cyclesPerAverage;
  int cycleNumber;
  double pa_q(vector<double> &f, vector<double> &cdf, int &j);
  vector<double> gvalArray;
  vector<double> cdfArray;
  bool checkConvergenceFlag;
  bool verifyConvergence;
  bool checkConvergence();
  double entropyConvergence;
  bool iterateGamma;

  bool levelSampling;
  bool unequalMasses;
  int activeBond; // bond that can form or break
  long int startingState;
  long int endingState;
  std::vector<bool> bondActive;
  std::vector<bool> stateActive;

  bool poolSample;
  int numPools;
  std::vector<ConfigurationPool> poolList; 
  ConfigurationPool globalPool;
  int numTraps;


  std::vector< std::vector<double> > reachedMin;
  bool recordAllDistances;

  bool dihedralOutput;
  bool histogramsInitialized;
  BiasHistograms* biasAngles;


  Matrix lab_frame(int index);
  void draw_atom(int index);

  void generateMolFile();
  void writeFrame();
  void outputMovie(int frames,double dt);
  void outputAcf();

  void initializePositions();
  void initializeVelocities();
  void initializeVelocitiesOld();
  int initPairs();

  void InitEventList(int num);
  void setupEvents();
  bool processEvent();
  bool processEventOld();
  double getNextEvent();
  double calculateCollisionTime(int i1,int i2, EventType& evType);
  void calculateAllCollisions(int na, int nb);
  void calculateCollisionPair(int i1, int i2);
  bool MDinitialized;

  void calculateAllBondCollisions();
  void calculateBondCollisions(int i1, int i2);
  void doNextCollision();
  void checkCollision(int collider1, int collider2, double cTime, double dsq);

  void resetSystem();
  void saveArrays();
  void restoreArrays();

  void ScheduleEvent(int idA, int idB, double tEvent, double deltaU, int collType, int bondIndex);
  void ScheduleEvent(Event& event);
  int insertInEventQ(int p);
  int insertInTree(int idNew);
  void processOverflowList();
  void deleteFromEventQ(int e);
  void DeleteEvent(int id);
  bool DrawNextEvent();
  void DrawNextEventQ();
  void replenishTree();
  void printNode(int id);
  void printList(int id);

  void printBondingState();

  bool debugMDflag;

  ofstream outFile, molFile, dihedralFile;
  void writeMovieVMD();
  bool generateImagesFlag,noAppend;
  std::vector<bool> wroteImage;
  std::vector<bool> wroteImage_save;
  
  void writeConfigurationVMD(int index);
  void generatePSF();
  void writeData();

  void findCrossings(int& num_crossings, int& writhe_i);

  void computeBareRates();
  void computeSinglePassageTime(int ib,bool output=false,bool writeToFile=true);
  void computePassageTime(int i);
  void computePassageTimeOld(int i);
  void testComputation(int i);
  void calcOuterFPT(double rc,long int starting, long int ending, std::vector<double> &data);
  void appendProblemState(ProblemPair &ps);


  void testFit(int ib);
  double computeLevelBareRates(SmoothPdfSpline &spdf);
  void plotPdf(int i, int j);
  double integrand(int i, int j, double x);
  double integrand(SmoothPdfSpline& spdf, double x);
  double trapzd(int i, int j, double a, double b, int n, double s);
  double trapzd(SmoothPdfSpline& spdf, double a, double b, int n, double s);
  double qsimp(int i, int j, double a, double b);
  double qsimp(SmoothPdfSpline& spdf, double a, double b);

  double qsimpOuter(SmoothPdfSpline& spdf, double a, double b);
  double qsimpInner(SmoothPdfSpline& spdf, double a, double b);
  double trapzdInner(SmoothPdfSpline& spdf, double a, double b, int n, double s);
  double trapzdOuter(SmoothPdfSpline& spdf, double a, double b, int n, double s);
  double integrandOuter(SmoothPdfSpline& spdf, double x);
  double integrandInner(SmoothPdfSpline& spdf, double x);


  double kramersRate(int i, int j);
  void printBias();

  std::vector<UniqueBond> uniqueBonds;
  void determineAllowedStates();
  double staircaseEntropy();
  double staircaseLevelEntropy(std::vector<int> &listActiveStates);
  double fptStaircase(int s1, int e1, double r1, double r2);
  double fptStaircase(int s1, int e1, double rc, double r1, double r2, std::vector<double> &data);
  double fptStaircaseBootstrap(int s1, int e1, double rc, double r1, double r2, std::vector<double> &data, int numBoot, double& stderr);
  void reduceEntropy();
  double reducedEntropyDifference(int sindex, int eindex);
  void writeEntropy();

  void runLevelStaircase(int sindex, int eindex, int bond,bool problemAnalysis=false);
  void samplePairStaircase();
  double staircaseLevelEntropy();

  bool updateEntropy(std::vector<int> & activeS);

  int findConnections(int i, int& aBond,double& fpt); // find best previous layer state to calculate entropy using fpt
  int findConnectionsVisits(int i, int& aBond); // using best previous layer using number of collisions
  void outputEntropy();
  void clearData();

  void testDihedrals(int startIndex);
  void biasGenerationMC();
  double redrawChain(std::vector<double> &bl, std::vector<double> &ba, std::vector<double> &phi, long int &configIndex);
  void drawAtom(Vector &r01_hat, Vector &r12_hat, double bl, double cos_alpha, double phi, Vector  &r23_hat,Vector& posAtom);
  double configBiasAtom(int biasIndex, int atomIndex, std::vector<double> &bl, std::vector<double> &cos_alpha, 
            std::vector<double>  &phiv, Vector &posAtom, double &new_ba, double &new_phi, bool select, bool& success);
  bool checkOverlap(int atomIndex, std::vector<Vector> &biasPos, Vector &posAtom, int& nb);
  bool checkOverlapEnergy(int atomIndex, std::vector<Vector> &biasPos, Vector &posAtom, double& pot);
  bool generateBiasMCTrial(int &biasStateIndex, std::vector<double> &bl, std::vector<double> &ba, std::vector<double> &phi,
                    double& logRatio,long int &trialState);
  void testInternals(std::vector<double> &new_ba, std::vector<double> &new_phi);
  double runBiasMC(int nsteps);
  int numBiasChoices;
  int numBiasStates; // number of histograms to load
  bool configBiasFlag;

  bool useWangLandau;
  int runTrajectoryWL(double time);
  void LevelWangLandau();
  void dumpEntropy();
  double relativeError(double& dev, double& err_i, int& errIndex);


  double extendedSamplePair(int start, int end, int aBond, double rmax,double &fpt, int activeLevels=0);
  void activatePreviousLevel(int state, std::vector<int> &previousLevelStates);
  bool extendedEntropy(std::vector<int> & map);
  double staircaseEntropyExtended();
  bool computeLevelEntropyNoFPT();

  void initializeFirstPool(int sIndex, ConfigurationPool &first);
  void readPoolFromFile(int index, ConfigurationPool &pool);
  void testPool(int intendedState, ConfigurationPool &sPool);
  void generateNextPool(int startingState,int endingState,ConfigurationPool &sPool, ConfigurationPool &ePool);
  void generatePoolFromPreviousLayers(int eState);
  void distributePool(ConfigurationPool &sPool, std::vector<int> &listActiveStates);
  void expandList(int listElement);
  void addDescendants(int list1, int list2);
  void runSingleState(ConfigurationPool &pool);
  
 public:


  ProteinChain(int nAtoms);
  ProteinChain(const string fileInputName, double betap);
  ~ProteinChain();

  void testCrankshaft();
  int crankshaftMC();
  bool annealMC();
  void runMC(int steps);
  void runOld(int nsteps,double outputInterval);
  void run();
  double computeEnergy(std::vector<Vector> &r, double& ke, double& pe);
  double computePotential(std::vector<double> &darray, long int& configNumber);
  double computePotential(std::vector<Vector> &r, long int& configNumber);
  double computeTrialPotential(std::vector<Vector> &r, long int& configNumber);
  double computeTrialPotential();

  void printOutput(double t);
  void checkpoint();
  void debugConfig();
  bool checkBonding();
  void PrintBondingDistances();
  void debugDynamics();

  long int configuration();
  long int configuration(std::vector<Vector> &trial);
  void loadConfiguration();
  void loadConfiguration(int index);
  void outputConfigurationCount();
  long int getConfigurationCount(int index) { return configurationCount[index]; }
  int runSampling(int numMC, double timeMD, double& mcRatio);
  void sample(int startIndex);
  void sampleLevels();
  void sampleSingleState(int index,bool = false);
  void saveStructure(int index,bool allowNewIndex = false);
  void levelSample(int startIndex,bool useWL = false);

  void whamEntropy();
  void computeEntropy();
  bool computeLevelEntropy(int& interval);
  bool computeLevelEntropyNew(int iteration);

  void initializeFPT();
  void computeLevelPassageTimes();

  double computeWritheAndGyration(double& gyr, double& maxD); // in analyze.cxx
 
  void scanFPT();
  void writeImage();
  void runModel();
  void runModelNew();
  void runStaircase(int index1, int eindex, int bond);
  void runSingleState();
  void runIsolatedStates();
  double analyzeDihedrals();
  void runProblemPairs();
  void runProblemPairsVariable();
  void runExtendedProblems();
  void runProblems();
  void WangLandau(int instance);
  void WangLandau(int instance,std::vector<int> &activeStates);
  int runWangLandauStates(int iterations, int iterationNumber, double& gamma, 
            std::vector<int> &states, std::vector<double> &Sval);

  double runTrajectory(int levelState1, int levelState2, double lambda, double time, int iterations=1);
  void generateState(int sIndex);

  void setBeta(double beta_p){ beta = beta_p; }
  double getBeta(){ return beta; }
};
