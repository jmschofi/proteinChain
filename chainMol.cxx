#include "chainMol.h"

//double INFINITY = 1.0E+12;
//double debugTime = 0.0;
//double debugTime = 2389420.0;

#define RESET_TIME
//#define TIMING
//#define PRINT_OVERLAP
//#define OLDQ

//#define PRINT_EVENT
//#define PRINT
//#define PRINT_BONDING
//#define PRINT_BONDING_CALC

//#define PRINT_BONDING_ENTRY
//#define PRINT_BONDING_LEAVE
//#define PRINT_MC

//#define PRINT_SW
//#define TEST_RESET
//#define TEST_RATES
//#define TEXTFILE
//#define PRINT_GVALUE
#define VERBOSE

ProteinChain::ProteinChain(const string inFileName, double betap){
  settings.read(inFileName);
  typedef std::chrono::high_resolution_clock myclock;
  myclock::time_point beginning = myclock::now();

  potEnergy = 0.0;  totEnergy = 0.0;
  potEnergy1 = 0.0; totEnergy1 = 0.0;
  kineticEnergy = 0.0;
  report.row("SYSTEM","");
 
  numAtoms = settings.getorset_int("numAtoms", 64);                report.COLUMN(numAtoms);
  numInteractions = 2*numAtoms - 3;
  pos.resize(numAtoms);// = new Vector[numAtoms];
  posTrial.resize(numAtoms);// = new Vector[numAtoms];
  posSave.resize(numAtoms);// = new Vector[numAtoms];
  vel.resize(numAtoms);// = new Vector[numAtoms];
  velSave.resize(numAtoms);// = new Vector[numAtoms];
  globalPool.setNumAtoms(numAtoms);

  isBondedAtom.resize(numAtoms);// = new bool[numAtoms];

  interactionNumber.resize(numAtoms);// = new int[numAtoms];
  localTime.resize(numAtoms);// = new double[numAtoms];
  timeNow = 0.0;
  globalTime = 0.0;

  for (int i=0;i<numAtoms;i++) {
    interactionNumber[i] = 0;
    localTime[i] = 0.0;
    isBondedAtom[i] = false;
  }
  //
  //  Initialization of the TREE
  poolSize = 500 * numAtoms ;
  evTree = NULL;
  scale = settings.getorset_int("scale" ,500);  report.COLUMN(scale); // used in hybrid queue to determine tree/linked list balance

#ifdef OLDQ
  evTree = new Event[poolSize+10];
  nlists = 500000;
  linearLists = new int[nlists+1];
  InitEventList(numAtoms);
#else
  hybridQ = new HybridQueue(scale); // test new hybrid queue
#endif
  
  mass = settings.getorset_DOUBLE("mass", 1.0);                        report.COLUMN(mass);
  
  //
  beta = 1.0; // not used
  double Etrans = 1.5/beta; report.COLUMN(Etrans); // average translational kinetic energy

  numAv = 0;
  Eav = 0.0;
  PEav = 0.0;
  KEav = 0.0;

  numIterations = settings.getorset_int("numIterations", 10); report.COLUMN(numIterations);
  numCycles = settings.getorset_int("numCycles",5); report.COLUMN(numCycles);
  cyclesPerAverage = settings.getorset_int("cyclesPerAverage",10); report.COLUMN(cyclesPerAverage);

  numMC = settings.getorset_int("numMC",100); report.COLUMN(numMC);
  configMCInterval = settings.getorset_int("configMCInterval", 300);
  report.COLUMN(configMCInterval);

  timeMD = settings.getorset_double("timeMD", 10000.0); report.COLUMN(timeMD);

  staircaseTime = settings.getorset_double("staircaseTime", 20.0); report.COLUMN(staircaseTime);
  isolationTime = settings.getorset_double("isolationTime",500000.0); report.COLUMN(isolationTime);

  configFlag = false;
  configFlag = settings.verify("configFlag"); report.COLUMN(configFlag);
  configInterval = settings.getorset_DOUBLE("configInterval", 0.1);
  if (configFlag) {
      maxNumDistances = settings.getorset_int("maxNumDistances", 2000000);
      report.COLUMN(configInterval);
      report.COLUMN(maxNumDistances);
  }
  verifyConvergence = settings.verify("verifyConvergence"); report.COLUMN(verifyConvergence);
  maxEntropyScaler = settings.getorset_DOUBLE("maxEntropyScaler", 0.1); report.COLUMN(maxEntropyScaler);
  minEntropyScaler = settings.getorset_DOUBLE("minEntropyScaler", 0.01); report.COLUMN(minEntropyScaler);
  entropyScaler = maxEntropyScaler;     

  poolSample = false; // only allow if doing level sampling

  iterateGamma = settings.verify("iterateGamma"); report.COLUMN(iterateGamma);

  double minBL = settings.getorset_DOUBLE("minBL", 1.0);                      report.COLUMN(minBL);
  double hardCore = settings.getorset_DOUBLE("hardCore", 1.25); report.COLUMN(hardCore);
  minBondLengthSq_hc = hardCore*hardCore; // hard core diameter
  innerBondDistance = hardCore;

  //minBondLengthSq_hc = minBL*minBL;
  minBondLengthSq1 = minBL*minBL;  // nearest neighbor bond min defines length unit
  double maxBL = settings.getorset_DOUBLE("maxBL", 1.2);                      report.COLUMN(maxBL);
  maxBondLengthSq1 = maxBL*maxBL;  // nearest neighbor bond max

  double minBL2 = settings.getorset_DOUBLE("minBL2", 1.7);                      report.COLUMN(minBL2);
  minBondLengthSq2 = minBL2*minBL2; 
  double maxBL2 = settings.getorset_DOUBLE("maxBL2", 2.0);                      report.COLUMN(maxBL2);
  maxBondLengthSq2 = maxBL2*maxBL2; // next nearest neighbor

  initialConfigFlag = settings.verify("initialConfigFlag");            report.COLUMN(initialConfigFlag);

  interactionFilename = settings.get_string("interactionFilename","bonds.dat"); report.COLUMN(interactionFilename);

  
  std::ifstream inputBondFile(interactionFilename);
  inputBondFile >> numBonds;
  report.COLUMN(numBonds);
  numBondsFixed = numBonds;
  maxBonds = 26;
  if (numBonds > numberBits){ 
    std::cerr << "  Max bonds for bitset analysis set to " << numberBits << " in compilation.  Must increase to at least " << numBonds << "." << std::endl;
    std::exit(1);
  }
 
  if (numBonds > 0){
    bondsInitialized = false;
    bondIndex1.resize(numBonds);// = new int[numBonds];
    bondIndex2.resize(numBonds);// = new int[numBonds];
    bondDistance.resize(numBonds);// = new double[numBonds];
    bondExists.resize(numBonds);// = new bool[numBonds];
    for (int i=0;i<numBonds;i++){
      inputBondFile >> bondIndex1[i] >> bondIndex2[i] >> bondDistance[i];
      isBondedAtom[ bondIndex1[i] ] = true;
	    isBondedAtom[ bondIndex2[i] ] = true; // these atoms in bonded interactions
	    bondExists[i] = false;
    }
 
  }


  dihedralOutput = settings.verify("dihedralOutput"); report.COLUMN(dihedralOutput);
  if (dihedralOutput)
  {
    string dihedralFilename = settings.get_string("dihedralFilename","dihedrals.dat"); report.COLUMN(dihedralFilename); 
    dihedralFile.open(dihedralFilename,std::ios::app);

  }

  histogramsInitialized = false;
  biasAngles = 0;

  pcrit = settings.getorset_DOUBLE("pcrit", 0.5);                      report.COLUMN(pcrit); 
  biasFlag = settings.verify("biasFlag"); report.COLUMN(biasFlag);
  if (biasFlag){

    biasFilename = settings.get_string("biasFilename","bias.dat"); 
    report.COLUMN(biasFilename);
    std::ifstream inputBiasFile(biasFilename);

	  numStates = 0;
	  std::string line;
	  while (std::getline(inputBiasFile, line)) ++numStates;
	  inputBiasFile.clear();
	  inputBiasFile.seekg(0, std::ios::beg);
    numStatesFixed = numStates;
    if (numStates != pow(2, numBonds)){
      std::cerr << "Number of states expected is " << pow(2,numBonds) << "  but biasing file " << biasFilename << " has data for "
	            	<< numStates << " states." << std::endl;
      std::exit(1);
    }
    
    Sbias.resize(numStates);//= new double[numStates];
    dS.resize(numStates);
    entropyMatrix.resize(numStates,std::vector<double>(cyclesPerAverage)); // 2d-array entropyMatrix[state][cycleIndex]

    std::vector<double> Smax(numBonds);// = new double[numBonds];
    for (int i=0;i<numBonds;i++) Smax[i] = -2.0;

    stateConnected.resize(numStates);
   
    for (int i=0;i<numStates;i++)
    {
        int state, decState;
        char stateLabel[126];
        double energy, entropy;
        inputBiasFile >> state >> energy >> entropy >> decState >> stateLabel;
        Sbias[decState] = entropy;
        stateConnected[i] = 0;
        //std::cout << "  decimal state " << decState << " had bias " << entropy << std::endl;
    }

   
    for (int i=1;i<numStates;i++) Sbias[i] -= Sbias[0]; // most bonded state has zero entropy;
    Sbias[0] = 0.0;   
    

    //for (int i=0;i<numBonds;i++){
     // std::cerr << "   Smax for " << i << " bonds is " << Smax[i] << std::endl;
   // }

    for (int i=0;i<numStates;i++){
      if (Sbias[i] < 0.0){
		    bitset<numberBits> integerIndex(i);
		    int eIndex = integerIndex.count();
		    //Sbias[i] = Smax[eIndexv
		    //std::cerr << "  Setting Sbias to " << Sbias[i] << "  for state " << i 
		    //	  << " with " << eIndex << " bonds." <<std::endl;
      }

    }
    //delete[] Smax;
   
  } else {
    //  no bias flag, first iteration
    numStates = pow(2, numBonds);
    numStatesFixed = numStates;
    Sbias.resize(numStates);// = new double[numStates];

    std::vector<double> Smax(numBonds);// = new double[numBonds];
    for (int i=0;i<numBonds;i++) Smax[i] = 3.0*(numBonds - i);
    entropyMatrix.resize(numStates,std::vector<double>(cyclesPerAverage)); // 2d-array entropyMatrix[state][cycleIndex]


    for (int i=0;i<numStates;i++){
      bitset<numberBits> integerIndex(i);
      int eIndex = integerIndex.count();
      double energy = -1.0*eIndex;
      Sbias[i] = 3.0*energy;
      
    }

    biasFlag = true;
  } // end of bias set up


 for (int i=0;i<numStates;i++) {
   State state_i(numAtoms, numBonds, i);
   states.push_back(state_i);
   states_save.push_back(state_i);
   states[i].setIndex(i);
   states_save[i].setIndex(i);
   for (int j=0;j<numBonds;j++) {
     states[i].setVisits(j,0);
     states_save[i].setVisits(j,0);
   }
  }

  problemStates.reserve(numStates);
  std::sort(states.begin(), states.end(), sortByLevel); // note that original indices for fpt data are stored in configuration
  std::sort(states_save.begin(), states_save.end(), sortByLevel);
  stateIndexToState.resize(numStates);
  stateIndexToState_save.resize(numStates);

  levelTarget.reserve(numStates);
  for (int i=0;i<numStates;i++) levelTarget[i] = false;
  if (dihedralOutput) levelTarget[numStates-1] = true; // record last states' dihedrals


  for (int i=0;i<numStates;i++){
        for (int s=0;s<numStates;s++) {
            if (states[s].getDecimalIndex() == i) {
              stateIndexToState[i] = s;
              stateIndexToState_save[i] = s;
            }
        }
  }

  //  print out mapping

  for (int i=0;i<numStates;i++)
  {
      states[i].setIndex(i);
      states_save[i].setIndex(i);
      //std::cout << "  ordered state index " << i << ": ";
      //states[i].print();
      //std::cout << " mapping for decimal index to ordered state index: " << i << " -> " << stateIndexToState[i] << std::endl;
  }


  allowedState.resize(numStates);// = new bool[numStates];  // used when staircase is active
  numAllowedStates = numStates;
  for (int i=0;i<numStates;i++) allowedState[i] = true;

  bondActive.resize(numBonds);// = new bool[numBonds];
  stateActive.resize(numStates);
  for (int i=0;i<numBonds;i++) bondActive[i] = true;//
  for (int i=0;i<numStates;i++) stateActive[i] = true; // default to all states active

  distanceMatrix.resize(numBonds);

  // obtain a seed from the timer
  myclock::duration d = myclock::now() - beginning;
  unsigned seed2 = d.count();
  generator.seed(seed2); /// initialize random generator

  recordAllDistances = false;
  recordAllDistances = settings.verify("recordAllDistances"); // attempt to find unreachable states
  if (recordAllDistances)
    {
        reachedMin.resize(numAtoms);
        for (int i=0;i<numAtoms;i++) 
        {
            reachedMin[i].resize(numAtoms);
            for (int j = 0; j < numAtoms;j++) reachedMin[i][j] = std::numeric_limits<double>::max();
        }
  }
 
  configBiasFlag = settings.verify("configBiasFlag");
  if (configBiasFlag)
  {
    numBiasChoices = settings.getorset_int("numBiasChoices", 10); report.COLUMN(numBiasChoices);
    numBiasStates = settings.getorset_int("numBiasStates",0); report.COLUMN(numBiasStates);
  }

  levelSampling = false;
  levelSampling = settings.verify("levelSampling"); report.COLUMN(levelSampling);
  if (levelSampling){
      activeBond = settings.getorset_int("activeBond",0);  report.COLUMN(activeBond);
      startingState = settings.getorset_int("startingState",0);  report.COLUMN(startingState); 
      bitset<numberBits> startingIndex(startingState);

      int eStart = startingIndex.count();

      startingIndex.flip ( activeBond );
      int eEnd = startingIndex.count();
    
      endingState = startingIndex.to_ulong();
      if (eStart > eEnd){
          // ending state has more bonds and lower entropy
          int temp = startingState;
          startingState = endingState;
          endingState = temp;
      }
      
      for (int i=0;i<numBonds;i++) bondActive[i] = false;
      bondActive[activeBond] = true;

      loadConfiguration(startingState);

      initialConfigFlag = false;


      std::bitset<numberBits> bitState1(startingState);
      std::bitset<numberBits> bitState2(endingState);
      std::vector<bool> state1;
      std::vector<bool> state2;
      for (int j=0;j<numBonds;j++) {
          state1.push_back( (bool)bitState1[j] );
          state2.push_back( (bool)bitState2[j] );
      }

      std::cout << "  Starting in state " << currentState <<  " states are " << startingState << " - " << endingState << std::endl;
      std::cout << "   Bonding patterns " << state1 << " " << state2 << std::endl;
      std::cout << "  Starting entropy values are " << Sbias[startingState] << " and " << Sbias[endingState] << std::endl;
      poolSample = settings.verify("poolSample"); report.COLUMN(poolSample); // only useful if level sampling
      if (poolSample)
      {
        numPools = settings.getorset_int("numPools", 1); report.COLUMN(numPools);
      }

  }

  zeroState = settings.getorset_int("zeroState",numStates-1);  report.COLUMN(zeroState); 
  randomSeed = settings.getorset_int("randomSeed",0); report.COLUMN(randomSeed);
  rndinit(0, randomSeed);

  if (levelSampling == false){ 
    endingState = numStates-1; // last state is target stateß
    if (initialConfigFlag){
      loadConfiguration();
    } else {
      initializePositions();
      initializeVelocities();
    }
  }

  redrawVelocityFlag = settings.verify("redrawVelocityFlag");            report.COLUMN(redrawVelocityFlag);
  if (redrawVelocityFlag) {
    redrawVelocityInterval = settings.getorset_DOUBLE("redrawVelocityInterval", 50.0); report.COLUMN(redrawVelocityInterval);
    redrawVelocityTime = redrawVelocityInterval;
    initializeVelocities();
  }

  spdf = 0;
  dataFlag = settings.verify("dataFlag");                            report.COLUMN(dataFlag);
  if (dataFlag){
    dataInterval = settings.getorset_DOUBLE("dataInterval", 50.0); report.COLUMN(dataInterval);
    dataTime = dataInterval;

  }

  useWangLandau = false;
  useWangLandau = settings.verify("useWangLandau");report.COLUMN(useWangLandau);

  movieFlag = settings.verify("movieFlag");                            report.COLUMN(movieFlag);

  if (movieFlag){
    movieInterval = settings.getorset_DOUBLE("movieInterval", 50.0); report.COLUMN(movieInterval);
    movieTime = movieInterval;
    movieInitialized = false;
    generateMolFile();
    outFile.open("chain.plb");
    float temp=numAtoms;
    outFile.write((char *)(&temp), sizeof(float));
    temp=0.0;
    outFile.write((char *)(&temp), sizeof(float));  
  }

  checkpointInterval = settings.getorset_DOUBLE("checkpointInterval",5);            report.COLUMN(checkpointInterval);
  outputInterval = settings.getorset_DOUBLE("outputInterval",0.1);     report.COLUMN(outputInterval);

  syncInterval = settings.getorset_DOUBLE("syncInterval",100.0);     report.COLUMN(syncInterval); 

  errorInterval = 1232.0;
  timeError = errorInterval; // TEST of reset function

  acfFlag = settings.verify("acfFlag");report.COLUMN(acfFlag);
  if (acfFlag) {
    acfInterval = settings.getorset_DOUBLE("acfInterval", 0.1);
    acfTotalInterval  = settings.getorset_DOUBLE("acfTotalInterval", 10.0);
    acfCorrelationTime = acfTotalInterval;
    acfPoints = 0;
    vij_rel_av = 0.0;
    acf = new Correlation<double,double>;
    acf->setup(int(acfTotalInterval/acfInterval), acfInterval, acfCorrelationTime); 
  }
  generateImagesFlag = false; 
  generateImagesFlag =  settings.verify("generateImagesFlag"); report.COLUMN(generateImagesFlag);
  wroteImage.resize(numStates);
  for (int i=0;i<numStates;i++) wroteImage[i] = false;

  if (generateImagesFlag==true){
    noAppend = false;
    noAppend = settings.verify("noAppend"); report.COLUMN(noAppend);
  }

  errorFlag = false;
  equilibrateFlag = settings.verify("equilibrateFlag"); report.COLUMN(equilibrateFlag);
  runTime = settings.getorset_DOUBLE("runTime",1.0);            report.COLUMN(runTime);
  basicTimeMD = 1.0;
  basicTimeMD = settings.getorset_DOUBLE("basicTimeMD",1.0); report.COLUMN(basicTimeMD);
  


  if (outputInterval > runTime) outputInterval = runTime;
  checkpointTime = checkpointInterval; report.COLUMN(checkpointTime);
  
  //   Start timing
  ret0=gettimeofday(&tv,&tz);
  sec0=(double )tv.tv_sec;
  microsec0=0.000001*( (double )tv.tv_usec);
  elapsedTime = 0.0;
  picoTime = 0.0;
  counterPicoTime = 0;

  numBondEvents = 0;
  numBondMC = 0;
 


  numStates = 1;
  if (numBonds > 0){
    numStates = pow(2, numBonds);
    configurationCount.resize(numStates);// = new long int[numStates];
    
    for (int i=0; i < numStates; i++) configurationCount[i] = 0;


    std::cout << std::endl << "#       " << numBonds << " Bonding interactions and " << numStates << " states." << std::endl;
    std::cout << "#   i1-i2\t\tdistance" << std::endl;
    std::cout << "# --------------------------" << std::endl;
    for (int i=0;i<numBonds;i++){
      std::cout << "#";
      if (bondActive[i]){
          std::cout << " T ";
      } else if (bondExists[i]) {
          std::cout << " F ";
      } else {
          std::cout << " O ";
      }
 
      std::cout << std::setw(3) << std::right << bondIndex1[i] << "-" << bondIndex2[i] << "\t\t" << bondDistance[i] << std::endl;
    }
  }
  std::cout << "#   ________________________________" <<  std::endl << std::endl;
  determineAllowedStates();

  collider1 = collider2 = -1;
  if (equilibrateFlag){
    double saveRunTime = runTime;
    bool saveDataFlag = dataFlag;
    bool saveConfigFlag = configFlag;
    configFlag = false;
    double saveVelocityInterval = redrawVelocityInterval;
    redrawVelocityInterval = 2000.0;
    redrawVelocityTime = redrawVelocityInterval;
    initializeVelocities();
    dataFlag = false;
    equilibrationTime = settings.getorset_DOUBLE("equilibrationTime",100.0);report.COLUMN(equilibrationTime);
    runTime = equilibrationTime;
    potEnergy = computePotential(pos, currentState);
    initialEnergy = computeEnergy(pos, kineticEnergy1,potEnergy1);
    totEnergy1 = initialEnergy;
    totEnergy = totEnergy1;
    potEnergy = potEnergy1;
    outputTime = outputInterval;
	
    std::cout << std::endl << "#   ________________ Equilibration Period ________________    output Time = " <<  outputTime <<  std::endl << std::endl;
    run(); // equilibration run
    // reset parameters

	
    initializeVelocities();
    timeNow = 0.0;
    globalTime = 0.0;
    runTime = saveRunTime;
    outputTime = outputInterval; 
    checkpointTime = checkpointInterval; 
    redrawVelocityInterval = saveVelocityInterval;
    redrawVelocityTime = redrawVelocityInterval;
    dataFlag = saveDataFlag;
    configFlag = saveConfigFlag;
    for (int i=0;i<numAtoms;i++) localTime[i] = timeNow;

    initialEnergy = computeEnergy(pos, kineticEnergy1,potEnergy1);
    totEnergy1 = initialEnergy;
    totEnergy = totEnergy1;
    potEnergy = potEnergy1;
    numAv = 0;
    Eav = 0.0;
    PEav = 0.0;
    KEav = 0.0;

    std::cout << "#                       Done Equilibration" << std::endl;
  }
  
  debugMDflag = settings.verify("debugMDflag");
  if (debugMDflag){
	  report.COLUMN(debugMDflag);
	  outputInterval = 1.1;
	  outputTime = outputInterval;
	  runTime = timeNow + timeMD;
	  std::cerr << "  Will debug MD for " << timeMD << " timesteps" << std::endl;
	  debugDynamics();
	  std::exit(1);
  }
  

  outputTime = outputInterval;
  setupEvents();
  std::cout << "#                      Beginning run of length " << runTime << std::endl;
  std::cout << "#      time \t E \t   Ke      PE  PE1 \t  Eav\t    Ke_av(/kT)\tPE_av " << endl;
  std::cout << "#   ______________________________________________________________________________" <<  std::endl;
  
  potEnergy = computePotential(pos, currentState);
  initialEnergy = computeEnergy(pos, kineticEnergy1,potEnergy1);
  totEnergy1 = initialEnergy;
  totEnergy = totEnergy1;
  potEnergy = potEnergy1;
  printOutput(timeNow);
}




ProteinChain::~ProteinChain(){

  //outputConfigurationCount();

  //if (dataFlag) computeBareRates();
  //if (dataFlag) computeLevelPassageTimes();

  //delete[] pos;
  //delete[] posTrial;
  //delete[] posSave;
  //delete[] vel;
  //delete[] velSave;
  //delete[] isBondedAtom;
  //delete[] interactionNumber;
  //delete[] localTime;
  //delete[] allowedState;

  //if (biasFlag) delete[] Sbias;
  if (acfFlag) delete acf;
  //if (generateImagesFlag) delete[] wroteImage;
  //if (numBonds > 0) {
	  //delete[] bondIndex1;
	  //delete[] bondIndex2;
	  //delete[] bondDistance;
	  //delete[] bondExists;
  //}

  if (spdf){
    for (int i=0;i<numStates;i++){
      delete[] spdf[i];
    }
    delete spdf;
  }

  //if (numBonds > 0) delete[] configurationCount;
  outFile.close();

  if (biasAngles != 0) delete biasAngles;

  if (hybridQ)  delete hybridQ;
}

void ProteinChain::determineAllowedStates(){

    //std::cerr << "  Determining the number of allowed states of the " << numStates << " bit states." << std::endl;
    std::vector< std::pair<int,int> > duplicatedPair; // index of bonds with more than one bond
    uniqueBonds.push_back( UniqueBond(bondIndex1[0], bondIndex2[0], bondDistance[0], 0) );

    for (int i=1;i<numBonds;i++){
      Bond b(bondIndex1[i], bondIndex2[i], bondDistance[i],i);
      int foundIndex = -1;
      for (int k=0;k<uniqueBonds.size(); k++) {
          if ( uniqueBonds[k].duplicateBond(b) ) foundIndex = k;
      }
      if (foundIndex > -1){
        uniqueBonds[foundIndex].addBond(b);
        //std::cout << "  adding bond " << b.index << " to unique bond " << foundIndex << std::endl;
      } else {
        uniqueBonds.push_back( UniqueBond(b) );
        //std::cout << "  Unique bond " << b.index << " started:" << std::endl;
      }

      for (int j=0;j<i;j++) {
        if ( (bondIndex1[i] == bondIndex1[j]) and (bondIndex2[i] == bondIndex2[j]) ) {

         if (bondDistance[i] > bondDistance[j]) 
           duplicatedPair.push_back( std::make_pair(i,j) ); // bond i must be set when bond j is set 
         else
           duplicatedPair.push_back( std::make_pair(j,i) ); // greater bond distance is redundant bond

        } 

     } 
    }
     //
    std::sort(duplicatedPair.begin(), duplicatedPair.end());
    //
  
    for (int i=0;i<numStates;i++) allowedState[i] = true;
    for (int k=0; k< duplicatedPair.size(); k++){

      int bond1 = duplicatedPair[k].first; // forbidden bond when second index is not set
      int bond2 = duplicatedPair[k].second;

      //std::cout << "   Pairs are " << bond1 << " - " << bond2 << std::endl;

      for (int i=0;i<numStates;i++){
        if (allowedState[i] == false) continue;
        
        std::bitset<numberBits> bitState(i);

	      if ( !bitState[bond1] and bitState[bond2] ) {
          allowedState[i] = false;
          //std::cout << "# Assessing that State " << i << " with bonding pattern " << bitState.to_string() << " is not allowed." << endl;
        } else {
          //std::cout << "# Assessing that State " << i << " with bonding pattern " << bitState.to_string() << " is still allowed." << endl;
        }


      }
    }
    numAllowedStates = 0;
    for (int i=0;i<numStates;i++){
      if (allowedState[i]) {
        std::bitset<numberBits> bitState(i);
        //std::cout << "  State " << i << " with bonding pattern " << bitState.to_string() << " allowed." << endl;
        numAllowedStates++;
      }
    }


   std::cout << "# num allowed states is " << numAllowedStates << std::endl;
}

void ProteinChain::reduceEntropy(){
    //
    //  Routine computes entropy of reduced model in which only shortest bond length for each bonding interatction of unique beads is used
    //      For example, if beads 2-8 have rc = 1.5 and rc = 4, only the bond defined as rc=1.5 is used to define a state
    //         The entropy for the new state must include contributions from 1.5 < r < 4 and r>4.


    int numBondsReduced = uniqueBonds.size();
    int numReduced = std::pow(2,numBondsReduced);

    std::vector<double> combinedEntropy(numReduced);
    combinedEntropy[numReduced-1] = 1.0; // entropy of last state is zero

    for (int i=numStates-2;i>=0;i--){
      std::bitset<numberBits> bitState(i); // original indexing
      std::bitset<numberBits> toState;
      if (!allowedState[i] ) continue;  // no contrib to any entropy
      for (int j=0;j<numBonds;j++){
        if (bitState[j]){
          // original bond j is on
          for (int k=0;k<uniqueBonds.size();k++){
            if (uniqueBonds[k].containsIndex(j) ) toState.flip(k);
          }
        }
      }
      int entropyIndex = toState.to_ulong();
      
      std::cout << "  State  " << i << " with pattern " << bitState.to_string() << " contributes to entropy " << entropyIndex << std::endl;
      combinedEntropy[entropyIndex] += exp( float(Sbias[i]) );

    }

    std::cout << "  Final entropies for reduced model:" << std::endl;
    for (int i=0;i<numReduced;i++){
      std::cout << " S[" << i << "] = " << log( float(combinedEntropy[i]) ) << std::endl;
    }

}

void ProteinChain::initializePositions(){
  double initLength = sqrt( float(minBondLengthSq1) )+0.50*( sqrt( float(maxBondLengthSq1) ) - sqrt( float(minBondLengthSq1) ) );
  pos[0] = Vector(0.0,0.0,0.0);
  pos[1] = Vector(initLength,0.0,0.0);
  double dx=initLength*cos(float(M_PI/2.0) );
  double dy = initLength*sin( float(M_PI/2.0) );
  pos[2] = Vector(initLength+dx,dy,0.0);
  if (numAtoms < 3) return;
  for (int i=3;i<numAtoms;i++) draw_atom(i);


  Vector R(0,0,0);
  for (int i=0;i<numAtoms;i++) R += pos[i];
  R /= double(numAtoms);
 
  for (int i=0;i<numAtoms;i++) {
      pos[i] -= R;
      //std::cout << "  Initial position of bead " << i << " is " << pos[i] << std::endl;
  }

}

void ProteinChain::draw_atom(int index){
  if (index < 3) return;
  int prev_index = index-1;
  Matrix T_lab=lab_frame(prev_index);
  double mag = sqrt(minBondLengthSq1)+0.50*( sqrt(maxBondLengthSq1) - sqrt(minBondLengthSq1) );
  double cos_alpha = cos(M_PI/2.0);
  double phi = M_PI;
  double sin_alpha=sqrt(1.0-cos_alpha*cos_alpha);
  Vector dl(-mag*cos_alpha,-mag*sin_alpha*cos(phi),mag*sin_alpha*sin(phi));
  Vector dv_i=T_lab*dl;

  pos[index] = pos[index-1] + dv_i;

#ifdef PRINT
  // std::cout << "  Position of atom " << index << " is " << pos[index] << std::endl;
#endif
  
  return;
}


Matrix ProteinChain::lab_frame(int index){
    int prev_index = index-1;
    int prev2_index = index-2;
    Vector U1 =  pos[index] - pos[index-1];
    Vector u1= U1/U1.nrm();

    Vector U0 = pos[index-1] - pos[index-2];
    Vector u0 = U0/U0.nrm();

    Vector U2 = u1^u0;
    Vector z = U2/U2.nrm();
    Vector y = z^u1;
    Matrix T1lab;
    T1lab.setColumn(0,u1);
    T1lab.setColumn(1,y);
    T1lab.setColumn(2,z);

    double norm=1.0/T1lab.det();
    Vector x1=u1*norm;
    Vector y1=y*norm;
    Vector z1=z*norm;

    Matrix result;
    result.setColumn(0,x1);
    result.setColumn(1,y1);
    result.setColumn(2,z1);
    return result;
}
void ProteinChain::initializeVelocities()
{
  double pxtot=0.0;
  double pytot=0.0;
  double pztot=0.0;
  double sigma = sqrt(1.0/mass/beta); // standard deviation
  std::normal_distribution<double> distribution(0.0, sigma);
  for (int i=0;i<numAtoms;i++)
  {
    double px = distribution(generator);
    double py = distribution(generator);
    double pz = distribution(generator);
    pxtot += px;
    pytot += py;
    pztot += pz;

    vel[i] = Vector(px,py,pz);
  }

  double dpx = pxtot/numAtoms;
  double dpy = pytot/numAtoms;
  double dpz = pztot/numAtoms;
  pxtot = 0.0;
  pytot = 0.0;
  pztot = 0.0;

  for (int j=0;j<numAtoms;j++)
  {
    double pxnew = vel[j].x - dpx;
    double pynew = vel[j].y - dpy;
    double pznew = vel[j].z - dpz;

    pxtot += pxnew;
    pytot += pynew;
    pztot += pznew;

    vel[j] = Vector(pxnew,pynew,pznew);
  }

  MDinitialized = false;

}
void ProteinChain::initializeVelocitiesOld(){
  //
  //         First, set the bead momenta
  //
  //
  double u1,u2,Rv,px,py,pz;
  double pxtot=0.0;
  double pytot=0.0;
  double pztot=0.0;

  for (int i=0 ;i < numAtoms;i++){
    double width_i = 0.50*beta/mass;
    u1=rnd();
    u2=rnd();
    Rv=sqrt( -log( 1.0-u1 )/width_i );
    px=Rv*cos(2.0*M_PI*u2);
   
    u1=rnd();
    u2=rnd();
    Rv=sqrt( -log( 1.0-u1 )/width_i );
    py=Rv*cos(2.0*M_PI*u2);
   
    u1=rnd();
    u2=rnd();
    Rv=sqrt( -log( 1.0-u1 )/width_i );
    pz=Rv*cos(2.0*M_PI*u2);
    
    pxtot += px;
    pytot += py;
    pztot += pz;

    vel[i] = Vector(px,py,pz);
  }

  double dpx = pxtot/numAtoms;
  double dpy = pytot/numAtoms;
  double dpz = pztot/numAtoms;

  pxtot = 0.0;
  pytot = 0.0;
  pztot = 0.0;

  double pxnew,pynew,pznew;

  for (int j=0;j<numAtoms;j++){
    pxnew = vel[j].x - dpx;
    pynew = vel[j].y - dpy;
    pznew = vel[j].z - dpz;

    pxtot += pxnew;
    pytot += pynew;
    pztot += pznew;

    vel[j] = Vector(pxnew,pynew,pznew);
  }

  //Vector V(0,0,0);
  //for (int i=0;i<numAtoms;i++) V += vel[i];

  //cout << "The total center of mass velocity is " << V  << endl;

  //cout << "  Velocities initialized with beta = " << beta << endl;
  MDinitialized = false;
}

void ProteinChain::sample(int startIter){

  double maxMDTime = 500.0;
  double minMDTime = 4.0;
  maxMDTime = settings.getorset_double("maxMDTime", 500.0); report.COLUMN(maxMDTime);
  minMDTime = settings.getorset_double("minMDTime", 5.0); report.COLUMN(minMDTime);
  //timeMD = settings.getorset_double("timeMD", 5.0); report.COLUMN(timeMD);
  if (timeMD < minMDTime) timeMD = minMDTime;
  configInterval = timeMD + 0.01; // only one config per propagation time: turned off

  int mcConfigs = numMC/configMCInterval;
  int mdConfigs = timeMD/configInterval;
  configsPerState = settings.getorset_int("configsPerState",400);  
  //report.COLUMN(configsPerState);

  
  
  int configsIteration = mdConfigs;

  // if (configsIteration == 0) {
  //     std::cerr << " No configurations generated per iteration.  Reduce configIntegral <= timeMD" << std::endl;
  //     std::exit(1);
  // }
  
  int updateInterval = configsPerState * numAllowedStates;
  if (levelSampling) updateInterval = configsPerState * 2.0;

  int updateConfiguration = updateInterval;
  

  // int outputInterval = 10;
  outputTime = outputInterval;
  int sIndex = stateIndexToState[startingState];
  
  for (int i=0;i<numStates;i++)
  {
    for (int b=0;b<numBonds;b++)
    {
      states[i].setVisits(b,0); // zero counters
    }
  }

  vector<int> energyV;
  for (int i=0;i<numStates;i++){
      bitset<numberBits> integerIndex(i);
      energyV.push_back( -integerIndex.count() );
      //std::cerr << "  State " << i << " has energy " << energyV[i] << std::endl;
  }

  int mdCount = 0;
  double mcRatio = 0.0;
  double mcRatioAv = 0.0;
  Stopwatch timer;
  timer.start();

  double totalTime = 0.0;
  int totalMC = 0;
  

  numIterations = numCycles*updateInterval;
  cycleNumber = 0;


  std::cout << " propagation time is " << timeMD << " configInterval = " << configInterval << " updateInterval = "
        << updateInterval << " update configuration = " << updateConfiguration 
        << " iter = " << startIter << " start = " << startingState << "  end = " << endingState 
        << " timeNow = " << timeNow << " numIterations = " << numIterations << std::endl;

  for (int i=0;i<numStates;i++) configurationCount[i] = 0;
  localHist.clear();
  configInterval = timeMD - 0.01;
  syncInterval = timeMD + 1.0; // should'nt need syncing since reset after md trajectory 
  doneScanIteration = false;
  int poolIndex = 0;
  
  int numTrajectories = int(timeMD/basicTimeMD + 1); // number of short-trajectories of time length basicMDTime per sample
  
  //for (int  iter = startIter; iter <= numIterations; iter++)
  int iter = 0;
  while (doneScanIteration == false)
  {

    iter++;
    if (poolSample)
    {
      globalPool.drawFromPool(pos,poolIndex, generator); // get new state, index poolIndex in list
      currentState = configuration(pos);
    }

	for (int traj = 0; traj < numTrajectories;traj++) mdCount += runSampling(numMC, basicTimeMD, mcRatio);
	localHist[currentState]++; // not an event, but done each iteration
	
  if (poolSample)
  {
    globalPool.replaceInPool(poolIndex, pos);
  }

  totalTime += numTrajectories * basicTimeMD;
  totalMC += numMC; 
  
  mcRatioAv += mcRatio;

  if (iter >= updateConfiguration){
    
    // std::cerr << "  Number of events in tree = " << numInTree << "  and number in lists = "
    // 			<< numInLists << "  total = " << numInTree + numInLists << "  # elements = " << localHist.size() << std::endl;
    

    for (auto p : localHist) configurationCount[p.first] += p.second; 
    localHist.clear();  // copy local hist into count then clear
    
    std::ofstream configFile;
    configFile.open("configurations.dat");
    int numRealStates = 0;

    int totalN = 0;
    for (int state=0; state < numStates; state++){ 
      if (configurationCount[state] > 0) numRealStates++;
      totalN += configurationCount[state];
    }
    double meanCount = double(totalN)/numRealStates;
    for (int state = 0; state < numStates; state++)  {
      if (configurationCount[state] > 0) {
        configFile << state << " " << configurationCount[state]/meanCount << " " << energyV[state]  << std::endl;
      }
    }
    configFile.close();

    if (levelSampling) {
      if (checkConvergenceFlag){
        checkConvergence();
      } else {
        if (fabs( Sbias[startingState] - Sbias[endingState]) > 8.0)
        {
          std::cout << "  Inefficient due to large entropy difference:  starting staircase." << std::endl;
          runStaircase(startingState,endingState, activeBond);
          doneScanIteration = true;
        }
        else
        {
          bool done = computeLevelEntropy(updateInterval);
        }      
      }
    } else {
          computeEntropy();
    }
    
    int nVisits = states[sIndex].getVisits(activeBond);

    

    double bondTime = 100.0; 
    if (numBondEvents > 0) bondTime = double(totalTime)/numBondEvents;
  

    timeMD = 2.0*bondTime + 1.0; 
    
    if (timeMD < minMDTime) timeMD = minMDTime;
    // if (nVisits < configsPerState/5) {
    //   timeMD += 50.0; // make longer when number visits is too small
    //   std::cout << " nVisits = " << nVisits << "  for state " << startingState << " so setting long md time." << std::endl;
    //  }  else if ( nVisits < configsPerState/3){
    //    timeMD += 10.0;
    //  }

    if (timeMD > maxMDTime) 
    {
      std::cout << " Setting timeMD to max value since bondTime is " << bondTime << std::endl;
      timeMD = maxMDTime;
    }

    configInterval = timeMD + 0.01; // turn off
    syncInterval = timeMD + 1.0; // should'nt need syncing since reset after md trajectory 
    numTrajectories = int(timeMD/basicTimeMD + 1);
    states[sIndex].setVisits(activeBond,0);

    std::cout << " -----    Output of configuration counts for num recorded states = "
        << numRealStates << " (" << double(numRealStates)/numStates * 100
        << "%)" << std::endl; 
#ifdef VERBOSE
      std::cout << " trajectories/sample = " << numTrajectories << "     md acceptance ratio = " 
        << double(mdCount)/double(iter+1) << "   mc ratio = "
              << mcRatioAv/double(iter+1) << std::endl << "      bond event rate = " << bondTime;
#endif

    if (numBondEvents > 0) std::cout << "   time to break or form bond = " << bondTime
                              << " md time is " << timeMD << "   ";
    if (totalMC > 0) std::cout  << "  bond mc event = " << double(numBondMC)/double(totalMC) << std::endl;


    numBondEvents = 0;
    totalTime = 0.0;

    //if (numBondEvents > 0 and configInterval < 1.1*bondTime) {
      //   configInterval = timeMD - 0.1;
      //   syncInterval = timeMD + 1.0; // should'nt need syncing since reset after md trajectory 

      //   std::cerr << "       Updated propagation time to " << timeMD << std::endl;
    // }

    updateConfiguration += updateInterval;
    checkpoint();
    timer.stop();
    
    timer.start();
	  }
  }
  std::cout << "   Completed " << numCycles << "   iter = " << numIterations << " without convergence." << std::endl;
  return;
}


int ProteinChain::runSampling(int numMC, double tMD, double& mcRatio){
	//
	//  Run trajectory and MC sampling from single starting configuration
	//

  int mdAccept = 1;
  bool done = false;
  errorFlag = false;
  errorCount = 0;
  mcRatio = 0;
  int accepts = 0;
  for (int i=0;i<numAtoms;i++) {
	  posTrial[i] = pos[i];
	  vel[i].zero();
  }

#ifdef PRINT_MC
  std::cout << "  Starting MC moves in state " << currentState << " num moves = " << numMC << std::endl;
#endif

  //mcRatio = runBiasMC(numMC); // configs not recorded in routine
  for (int i=0;i<numMC;i++){
    accepts +=  crankshaftMC();
  }

  if (numMC > 0) mcRatio = double(accepts)/double(numMC);
#ifdef PRINT_MC
  std::cout << "  Crankshaft acceptance ratio is " << mcRatio << std::endl;
#endif

  beta = 1.0;
  
#ifdef RESET_TIME
  for (int i=0;i<numAtoms;i++){
    double dt = timeNow - localTime[i];
    pos[i] += vel[i]*dt;
    localTime[i] = 0.0;
    interactionNumber[i] = 0; // NEW
  }

  globalTime += timeNow;
  timeNow = 0.0; // avoid round-off errors?
#endif

  initializeVelocities();

  double ke, pe;
  double Etest = computeEnergy(pos, ke, pe);
  double Sbefore = ke + Sbias[currentState];

  //checkpointInterval = timeMD+1.0; // disable checkpointing
  checkpointTime = timeNow + checkpointInterval;
  if (dataFlag) dataTime = timeNow + dataInterval; 
  equilibrateFlag = false;
  movieFlag = false;
  acfFlag = false;
  //outputTime = timeMD + 1.0; // no output or check of energy conservation
  outputTime = checkpointInterval;
 
  runTime = timeNow + tMD;
  //if (timeMD <= configInterval){
       //std::cerr << "  No samples taken since md time is " << timeMD << " <= " << configInterval << std::endl;
       //configInterval = timeMD - 0.000000001;
       //std::cerr << " Setting configuration sampling time to be slightly less than md time." << std::endl;
  //}

  bondsInitialized = false;
  MDinitialized = false;

  bitset<numberBits> stateIndex(currentState);
  configFlag = true;
  setupEvents(); // 

  for (int i=0;i<numAtoms;i++){
      velSave[i] = vel[i];
      posSave[i] = pos[i] + vel[i]*(timeNow-localTime[i]);
      saveState = currentState;
      startTime = timeNow;
  }
 
  
  while (!done) done = processEvent();

	  
  if (!checkBonding() ) errorFlag = true;

  int configNumber = configuration();
  if (currentState != configNumber){
	  std::cerr << "  In dynamical update, configuration from direct calc. = "
				<< configNumber << "  and dynamically set to " << currentState << std::endl;
	  errorFlag = true;
  }

  Etest = computeEnergy(pos, ke, pe);
  double Safter = ke + Sbias[currentState];

  //bitset<numberBits> endIndex(currentState);
  //std::cout << " ... now in state " << currentState << " (" << endIndex << " )" << std::endl;

  if ( fabs(Safter-Sbefore) > 1.0e-7) {
	  std::cerr << "Possible error in total Hamiltonian conservation since Safter = " << Safter
				<< "   Sbefore = " << Sbefore << "    ds = " << Safter - Sbefore << std::endl;
	  errorFlag = true;
  }


  if (errorFlag){
	  errorCount = 0;
	  //debugConfig();
	  errorFlag = false;
	  for (int i=0;i<numAtoms;i++) {
		  pos[i] = posSave[i]; // reset state to try again
		  timeNow = startTime;
		  localTime[i] = timeNow;
	  }
    loadConfiguration(startingState); 
	  currentState = startingState;
	  bondsInitialized = false;
	  std::cout << "  potential energy is ";
	  potEnergy = computePotential(pos, currentState);
	  std::cout << "  pe = " << potEnergy << "  state is " << currentState << std::endl;
          mdAccept = 0;

  } 

  
  return mdAccept;
}

void ProteinChain::run(){
  bool done = false;
  if (MDinitialized == false) setupEvents();
  while (!done) done = processEvent();
}


void ProteinChain::outputMovie(int frames,double dt){
  for (int i=0;i<frames;i++){
    double time=double(i)*dt;
    for (int j=0;j<numAtoms;j++){
      Vector p1=pos[j]+vel[j]*time;
      float x1=p1.x;
      float y1=p1.y;
      float z1=p1.z;

      outFile.write((char *)(&x1), sizeof(float));
      outFile.write((char *)(&y1), sizeof(float));
      outFile.write((char *)(&z1), sizeof(float));
    }
  }
  outFile.flush();
}

double ProteinChain::computeEnergy(std::vector<Vector> &r, double& kEnergy, double& pEnergy){

  bitset<numberBits> integerIndex;
  integerIndex.reset(); // sets all bits to zero

  pEnergy = 0.0;
  kEnergy = 0.0;
  if (bondsInitialized == false){
    for (int i=0;i<numBonds;i++) bondExists[i] = false;
  }

  if (!bondsInitialized){
      for (int i=0;i<numBonds;i++) bondExists[i] = false;
  }

  for (int i=0;i<numBonds;i++)
  {
      Vector rij = r[ bondIndex1[i] ] - r[ bondIndex2[i] ];
      double rijval = rij.nrm();
      if ( rijval < bondDistance[i] ) 
      {
	      if (!bondsInitialized) bondExists[i] = true;
        integerIndex.set(i,true);
        pEnergy -= 1.0;
      }
  }
  
  for (int i = 0; i < numAtoms;i++){
    kEnergy += 0.5*mass*vel[i].nrm2();

    Vector ri = r[i];
    for (int j = i+1;j < numAtoms; j++){
      Vector rij = r[j]-ri;
      double dijsq = rij.nrm2();

      if (dijsq < minBondLengthSq_hc and (j > i+2) ){
		       std::cerr << " computeEnergy(Vector* r, double& kEnergy, double& pEnergy) :Error in hard core interaction "
		        	<< i << "-" << j << "  dijsq = "
		        	<< dijsq << " is less than " << minBondLengthSq_hc << std::endl;
		        pEnergy = std::numeric_limits<double>::max();
      }

      if ( (j == i+1) and (dijsq < minBondLengthSq1) ){
	         std::cerr << " computeEnergy(Vector* r, double& kEnergy, double& pEnergy) :Error in bonded interaction "
		          	<< i << "-" << j << "  dijsq = " << dijsq << " is less than " << minBondLengthSq1 << std::endl;
	          pEnergy = std::numeric_limits<double>::max();
      }

     if ( (j == i+2) and (dijsq < minBondLengthSq2) ){
	       std::cerr << "computeEnergy(Vector* r, double& kEnergy, double& pEnergy): Error in bond angle interaction "
	        	   << i << "-" << j << "  dijsq = " << dijsq << " is less than " << minBondLengthSq2
	        	   << " diff = " << dijsq - minBondLengthSq2 << std::endl;
	       pEnergy = std::numeric_limits<double>::max();
      }

      if ( (j == i+1) and (dijsq > maxBondLengthSq1) ){
	         std::cerr << "computeEnergy(Vector* r, double& kEnergy, double& pEnergy): Error in bonded interaction "
	            	<< i << "-" << j << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq1
	            	<< " diff = " << dijsq - maxBondLengthSq1 << std::endl;
        	  pEnergy = std::numeric_limits<double>::max();
      }

      if ( (j == i+2) and (dijsq > maxBondLengthSq2) ){
	         std::cerr << "computeEnergy(Vector* r, double& kEnergy, double& pEnergy): Error in bond angle interaction "
              		<< i << "-" << j << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq2 << std::endl;
	          pEnergy = std::numeric_limits<double>::max();
      }


    } 

  }



  bondsInitialized = true;

  long int configNumber = integerIndex.to_ulong();
  pEnergy = Sbias[configNumber];

  double tEnergy = kEnergy + pEnergy;
  return tEnergy;
}





double ProteinChain::computePotential(std::vector<Vector> &r, long int& configNumber){
  bitset<numberBits> integerIndex;
  integerIndex.reset(); // sets all bits to zero
  double pEnergy = 0.0;
  for (int i = 0; i < numAtoms;i++){
    Vector ri = r[i];
    for (int j = i+1;j < numAtoms; j++){
      Vector rij = r[j]-ri;
      double dijsq = rij.nrm2();

      if (dijsq < minBondLengthSq_hc and (j > i+2) ){
		      std::cerr << "Overlap hard core interaction " << i << "-" << j
					        << "  dijsq = " << dijsq << " is less than " << minBondLengthSq_hc << std::endl;
		      return std::numeric_limits<double>::max();
      }

      if ( (j == i+1) and (dijsq < minBondLengthSq1) ){
		        std::cerr << "Overlap bonded interaction " << i << "-" << j
			        		<< "  dijsq = " << dijsq << " is less than " << minBondLengthSq1 << std::endl;
		        return std::numeric_limits<double>::max();
      }

     if ( (j == i+2) and (dijsq < minBondLengthSq2) ){
        std::cerr << "Overlap bond angle interaction " << i << "-" << j
				    << "  dijsq = " << dijsq << " is less than " << minBondLengthSq2 << std::endl;
		    return std::numeric_limits<double>::max();
      }

      if ( (j == i+1) and (dijsq > maxBondLengthSq1) ){
		      std::cerr << "Overlap bonded interaction " << i << "-" << j
			    		<< "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq1 << std::endl;
		      return std::numeric_limits<double>::max();
      }

      if ( (j == i+2) and (dijsq > maxBondLengthSq2) ){
		      std::cerr << "Overlap bond angle interaction " << i << "-" << j
				        	<< "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq2 << std::endl;
		     return std::numeric_limits<double>::max();
      }


    } 

  }

  if (!bondsInitialized) {
      for (int i=0;i<numBonds;i++) bondExists[i] = false;
  }

  for (int i=0;i<numBonds;i++){
      Vector rij = r[ bondIndex1[i] ] - r[ bondIndex2[i] ];
      double rv = rij.nrm2();
      if ( rv < bondDistance[i]*bondDistance[i] ) {
		  if (!bondsInitialized){
			  bondExists[i] = true;
		  }
		  pEnergy -= 1.0;
		  integerIndex.set(i,true);
      } 
  }
  bondsInitialized = true;
  configNumber = integerIndex.to_ulong();

  pEnergy = Sbias[configNumber]; // entropy used as potential
  
  //std::cout << "  configuration " << integerIndex << " is configuration number " << configNumber << std::endl;


  return pEnergy;
}

bool ProteinChain::checkBonding(){
    // check bonding state with dynamic state
    bool bondsOk = true;
    for (int i=0;i<numBonds;i++){
        int i1 = bondIndex1[i];
        int i2 = bondIndex2[i];
        double dt1 = timeNow - localTime[i1];
        double dt2 = timeNow - localTime[i2];
        Vector rij = pos[i1] + vel[i1]*dt1 - pos[i2] - vel[i2]*dt2;
        double rsq = rij.nrm2();
        if (rsq < bondDistance[i]*bondDistance[i]){
            if (!bondExists[i]){
                bondsOk = false;
                std::cerr << "  Error in dynamic bond accounting at bond " << i << " between " << i1 << "-" << i2
		                    << " where bond should exist since rsq = " << rsq << std::endl;
			exit(1);
	          }
	      } else {
	          if (bondExists[i]){
	              bondsOk = false;
		            std::cerr << "  Error in dynamic bond accounting at bond " << i << " between " << i1 << "-" << i2
		               << " where bond should not exist since rsq = " << rsq << std::endl;
			exit(1);
	          }
	      }
    }

    return bondsOk;

}

double ProteinChain::computeTrialPotential(std::vector<Vector> &r, long int& configNumber){
  bitset<numberBits> integerIndex;
  integerIndex.reset(); // sets all bits to zero
  double pEnergy = 0.0;
  for (int i = 0; i < numAtoms;i++){
    Vector ri = r[i];
    for (int j = i+1;j < numAtoms; j++){
      Vector rij = r[j]-ri;
      double dijsq = rij.nrm2();

      if (dijsq < minBondLengthSq_hc and (j > i+2) ){
#ifdef PRINT_OVERLAP
        std::cout << " Overlap for pair " << i << "-" << j << " dij = " << sqrt(dijsq) << " less than "
              << sqrt(minBondLengthSq_hc) << std::endl;
#endif
		    return std::numeric_limits<double>::max();
      }

      if ( (j == i+1) and (dijsq < minBondLengthSq1) ){
		    return std::numeric_limits<double>::max();
      }

     if ( (j == i+2) and (dijsq < minBondLengthSq2) ){
		    return std::numeric_limits<double>::max();
      }

	   if ( (j == i+1) and (dijsq > maxBondLengthSq1) ){
		   return std::numeric_limits<double>::max();
	   }

	   if ( (j == i+2) and (dijsq > maxBondLengthSq2) ){ 
	  	 return std::numeric_limits<double>::max();
	   }
    } 

  }

  if (!bondsInitialized){
      for (int i=0;i<numBonds;i++) bondExists[i] = false;

  }

  for (int i=0;i<numBonds;i++){
      Vector rij = r[ bondIndex1[i] ] - r[ bondIndex2[i] ];
      double rv = rij.nrm2();
      if ( rv < bondDistance[i]*bondDistance[i] ) {
		    if (!bondsInitialized)
        {
          bondExists[i] = true;
		    }
		   pEnergy -= 1.0;
       integerIndex.set(i,true);
       //std::cout << "  Bond " << bondIndex1[i] << "-" << bondIndex2[i] << " in trial is set." << std::endl;
      } 
  }
  bondsInitialized = true;
  configNumber = integerIndex.to_ulong();

  //std::cout << "  configuration " << integerIndex << " is configuration number " << configNumber << std::endl;

  if (biasFlag) pEnergy = Sbias[configNumber];
  return pEnergy;
}


double ProteinChain::computeTrialPotential(){

  std::cout << "  ----------  computeTrialPotential() called " << std::endl;

	double pEnergy = 0.0;
	for (int i = 0; i < numAtoms;i++){
		Vector ri = posTrial[i];
		for (int j = i+1;j < numAtoms; j++){
			Vector rij = posTrial[j]-ri;
			double dijsq = rij.nrm2();
			
			if (dijsq < minBondLengthSq_hc and (j > i+2) ){
				std::cout << "Overlap hard core interaction " << i << "-" << j << "  dijsq = "
				        << dijsq << " is less than " << minBondLengthSq_hc << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+1) and (dijsq < minBondLengthSq1) ){
				std::cout << "Overlap interaction " << i << "-" << j
				          << "  dijsq = " << dijsq << " is less than " << minBondLengthSq1 << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+2) and (dijsq < minBondLengthSq2) ){
				std::cout << "Overlap ba interaction " << i << "-" << j << "  dijsq = " << dijsq << " is less than " << minBondLengthSq2 << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+1) and (dijsq > maxBondLengthSq1) ){
				std::cout << "Overlap interaction " << i << "-" << j << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq1 << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+2) and (dijsq > maxBondLengthSq2) ){
				std::cout << "Overlap interaction " << i << "-" << j << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq2 << std::endl;
				return std::numeric_limits<double>::max();
			}
		} 
	}

 

	for (int i=0;i<numBonds;i++){
		Vector rij = posTrial[ bondIndex1[i] ] - posTrial[ bondIndex2[i] ];
		double rv = rij.nrm2();
		if ( rv < bondDistance[i]*bondDistance[i] ) {	
      //std::cout << "  Bond " << bondIndex1[i] << "-" << bondIndex2[i] << " is on in trial." << std::endl;
			pEnergy -= 1.0;
		}
	}

	return pEnergy;
}


double ProteinChain::computePotential(std::vector<double> &rd, long int& configNumber){
	bitset<numberBits> integerIndex;
	integerIndex.reset(); // sets all bits to zero
	double pEnergy = 0.0;
	for (int i = 0; i < numAtoms;i++){
		Vector ri = pos[i];
		for (int j = i+1;j < numAtoms; j++){
			Vector rij = pos[j]-ri;
			double dijsq = rij.nrm2();

			if (dijsq < minBondLengthSq_hc and (j > i+2) ){
				std::cout << "Overlap hard core interaction " << i << "-" << j
						  << "  dijsq = " << dijsq << " is less than " << minBondLengthSq_hc << " time = " << timeNow <<  std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+1) and (dijsq < minBondLengthSq1) ){
				std::cout << "Overlap bonded interaction " << i << "-" << j << "  dijsq = "
						  << dijsq << " is less than " << minBondLengthSq1 
						  << " diff = " << dijsq - minBondLengthSq1 << " time = " << timeNow << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+2) and (dijsq < minBondLengthSq2) ){
				std::cout << "Overlap bond angle interaction " << i << "-" << j
						  << "  dijsq = " << dijsq << " is less than " << minBondLengthSq2 << " time = " << timeNow << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+1) and (dijsq > maxBondLengthSq1) ){
				std::cout << "Overlap bonded interaction " << i << "-" << j
						  << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq1
						  << " time = " << timeNow << std::endl;
				return std::numeric_limits<double>::max();
			}

			if ( (j == i+2) and (dijsq > maxBondLengthSq2) ){
				std::cout << "Overlap bond angle interaction " << i << "-" << j
						  << "  dijsq = " << dijsq << " is greater than " << maxBondLengthSq2
						  << " time = " << timeNow << std::endl;
				return std::numeric_limits<double>::max();
			}


		} 
		
	}

	for (int i=0;i<numBonds;i++){
		Vector rij = pos[ bondIndex1[i] ] - pos[ bondIndex2[i] ];
		double rv = rij.nrm2();
		rd[i] = rv;
		if ( rv < bondDistance[i]*bondDistance[i] ) {
			pEnergy -= 1.0;
			integerIndex.set(i,true);
		}
	}

	configNumber = integerIndex.to_ulong();
  pEnergy = Sbias[configNumber];

  std::cout << "  configuration " << integerIndex << " is configuration number " << configNumber << std::endl;

	return pEnergy;
}

void ProteinChain::printOutput(double t){

  std::cout << std::fixed << " " << setw(10) << t << "\t"
			<< std::setprecision(5) << totEnergy1/double(numAtoms-1) << "\t " << std::setprecision(5)  
			<< kineticEnergy1/double(numAtoms-1);
  std::cout.unsetf ( std::ios::fixed );
  std::cout << std::setw(3)   << " " << potEnergy1 << " " << std::setw(3)  << potEnergy << setw(6) << " " << currentState
    << " " << previousState;
  if (numAv > 0){
    double norm = double((numAtoms-1)*numAv);
    std::cout << "\t" << setw(6) << Eav/norm << "\t     "
			  << setprecision(5) << KEav/norm*beta << "\t" << setw(4) << PEav/double(numAv);
  }    
  std::cout << " " << distances.size() << std::endl;

  //Vector R(0,0,0);
  //Vector V(0,0,0);
  //for (int i=0;i<numAtoms;i++) {
  //   R += pos[i];
  //   V += vel[i];
  //}
  //R /= double(numAtoms);
  //V /= double(numAtoms);
  //std::cout << "  center of mass position is R = " << R << "         V = " << V << std::endl;
 
  // if (potEnergy != potEnergy1) {
  //   std::cout << "Error at time " << timeNow << " in potential energy." << std::endl;
  //   long int config = configuration();
  //   resetSystem();
  //   //exit(1);
  // }

  // if ( fabs(totEnergy1-initialEnergy) > 0.0000001){
  //   std::cout << "Error at time " << timeNow << " in total energy " << totEnergy1 << " not equal to initial energy " << initialEnergy << std::endl;
  //   long int config = configuration();
  //   resetSystem();
  //   //exit(1);
  // }

  // simulate an error
#ifdef TEST_RESET
  if (timeNow > timeError) {
    timeError += errorInterval;
    std::cerr << "  ##### Simulating an error at time " << timeNow << std::endl;
    resetSystem();
  }
#endif
}

void ProteinChain::printBondingState(){
     std::cout << "  Bonding flags at time " << timeNow << std::endl;
     for (int i=0;i<numBonds;i++){
         std::cout << "  Bond " << i << " between " << bondIndex1[i] << "-" << bondIndex2[i] << " has state " << bondExists[i]
	           << std::endl;
     }
}

void ProteinChain::debugDynamics(){

  binfile data_file;
  std::cerr << "  Loading configuration to debug.." << std::endl;

  stringstream filenamebase;
  filenamebase << "debugConfig.bin";
  std::string tmp = filenamebase.str();
  const char* name = tmp.c_str();

  data_file.open(name, "r");


  double rx, ry, rz, vx, vy, vz;
  for (int j=0;j<numAtoms;j++){
    //
    //    Read in center of mass coordinates
    //
    data_file >> rx >> ry >> rz >> vx >> vy >> vz;
    Vector r(rx, ry, rz);
    Vector v(vx, vy, vz);
    pos[j] = r;
    vel[j] = v;
    localTime[j] = timeNow;
  }
  bondsInitialized = false;
  potEnergy = computePotential(pos, currentState);
  totEnergy1  = computeEnergy(pos, kineticEnergy1, potEnergy1);

  std::cout << "  Before call to configuration, state is " << currentState << "  e = " << potEnergy << std::endl;
  printBondingState();


  currentState = configuration();
  previousState = currentState;
  setupEvents();
  runTime = timeMD;
  ScheduleEvent(-1, -1, runTime, 0.0,  STOP, -1);
  hybridQ->ScheduleEvent(-1,-1, runTime, 0.0, STOP, -1, -1,-1);

  std::cout << "  Before dynamics and after event set up, state is " << currentState << std::endl;
  bool done = false;
  while (!done) done = processEvent();

}

void ProteinChain::debugConfig(){
  binfile data_file;

  stringstream filenamebase;
  filenamebase << "debugConfig.bin";
  std::string tmp = filenamebase.str();
  const char* name = tmp.c_str();

  data_file.open(name, "w");
  std::cout << "#******* writing debug configuration: timeNow = " << timeNow << "  saved config. at time "
            << startTime << "  time of problem =" << timeNow - startTime
			<< " state is " << saveState << "  iteration is " << saveIter << std::endl;
  for (int j=0;j<numAtoms;j++){
    //
    //    write out phase space coordinates
    //
    data_file << posSave[j].x << posSave[j].y << posSave[j].z;
    data_file << velSave[j].x << velSave[j].y << velSave[j].z;
  }
 

}

void ProteinChain::checkpoint(){

 if (generateImagesFlag){
    if (!wroteImage[currentState]) {
      saveStructure(currentState);
      writeConfigurationVMD(currentState);
      if (noAppend) wroteImage[currentState]= true;
      std::cout << "  Wrote structure " << currentState << std::endl;
    }
  }


  //   zero the center of mass

  Vector R(0,0,0);
  double totalM = 0.0;
  for (int i=0;i<numAtoms;i++){
      R += pos[i] * mass;
      totalM += mass;
  }

  for (int i=0;i<numAtoms;i++) pos[i] -= R;



  binfile data_file;
  data_file.open("config.bin", "w");
  //std::cout << "#******* checkpoint: ";
  for (int j=0;j<numAtoms;j++){
    //
    //    write out phase space coordinates
    //
    data_file << pos[j].x << pos[j].y << pos[j].z;
    data_file << vel[j].x << vel[j].y << vel[j].z;
  }

  data_file.flush();
  data_file.close();
 
  ret1=gettimeofday(&tv,&tz);
  sec1=(double )tv.tv_sec;
  microsec1=0.000001*( (double )tv.tv_usec);
  elapsedTime = sec1 - sec0 + microsec1 - microsec0;

  //std::cout << " time elapsed " << elapsedTime << "   potEnergy = " << potEnergy1 << "     beta = " << beta << std::endl;
}


void ProteinChain::loadConfiguration(){
  binfile data_file;
  data_file.open("config.bin", "r");
  double rx, ry, rz, vx, vy, vz;
  for (int j=0;j<numAtoms;j++){
    //
    //    Read in center of mass coordinates
    //
    data_file >> rx >> ry >> rz >> vx >> vy >> vz;
    Vector r(rx, ry, rz);
    Vector v(vx, vy, vz);
    pos[j] = r;
    vel[j] = v;
    localTime[j] = timeNow;
  }
  potEnergy = computePotential(pos, currentState);
  totEnergy1  = computeEnergy(pos, kineticEnergy1, potEnergy1);

}

void ProteinChain::PrintBondingDistances(){
  for (int i=0;i<numBonds;i++){
        int i1 = bondIndex1[i];
        int i2 = bondIndex2[i];
        double dt1 = timeNow - localTime[i1];
        double dt2 = timeNow - localTime[i2];
        Vector rij = pos[i1] + vel[i1]*dt1 - pos[i2] - vel[i2]*dt2;
        double rsq = rij.nrm2();
        double r = sqrt(rsq);
        std::cout << " Bond " << i1 << "-" << i2 << " is at distance " << r 
                  << " rc = " << bondDistance[i] << " bondExists = " << bondExists[i] << std::endl;
      
  }


}

void ProteinChain::loadConfiguration(int index){
  binfile data_file;


  stringstream filenamebase;
  filenamebase << "structure." << index << ".bin";
  std::string tmp = filenamebase.str();
  const char* name = tmp.c_str();

  data_file.open(name, "r");
  std::cout << "#***  Reading in configuration " << name << std::endl;

  double rx, ry, rz, vx, vy, vz;
  for (int j=0;j<numAtoms;j++){
    //
    //    Read in center of mass coordinates
    //
    data_file >> rx >> ry >> rz >> vx >> vy >> vz;
    Vector r(rx, ry, rz);
    Vector v(vx, vy, vz);
    pos[j] = r;
    vel[j] = v;
    localTime[j] = 0.0;
  }
  bondsInitialized = false;
  currentState = index;
  previousState = currentState;
  potEnergy = computePotential(pos, currentState);
  totEnergy1  = computeEnergy(pos, kineticEnergy1, potEnergy1);
  std::cout << "##  Loaded configuration " << index << "  =  " << currentState 
        << "  potEnergy = " << potEnergy << std::endl;
  if (index != currentState) {
    std::cerr << "  Warning: File contains configuration " << currentState << " rather than " << index << "  with " << numBonds << " bonds." << std::endl;
    PrintBondingDistances();
    exit(1);
  }
 }


void ProteinChain::saveStructure(int index,bool allowNewIndex){



  for (int i=0; i < numAtoms; i++){
      double dt = timeNow - localTime[i];
      posTrial[i] = pos[i] + vel[i]*dt;
  }


  //   zero the center of mass

  Vector R(0,0,0);
  double totalM = 0.0;
  for (int i=0;i<numAtoms;i++){
      R += posTrial[i] * mass;
      totalM += mass;
  }

  for (int i=0;i<numAtoms;i++) posTrial[i] -= R;

  long int configNumber = -1;
  double pot = computePotential(posTrial, configNumber);
  if ((allowNewIndex == false) and configNumber != index){
    std::cerr << "  Error in savestructure routine since config is " << configNumber << " and not " << index << std::endl;
    exit(1);
  }

  std::cout << "  ******************  Wrote configuration " << index << std::endl;

  binfile data_file;
  std::ofstream textFile;

  stringstream filenamebase;
  filenamebase << "structure." << index << ".bin";
  std::string tmp = filenamebase.str();
  const char* name = tmp.c_str();

  data_file.open(name, "w");

#ifdef TEXTFILE
  stringstream base2;
  base2 << "structure." << index << ".txt";
  std::string tmp2 = base2.str();
  const char* nametxt = tmp2.c_str();
  textFile.open(nametxt,std::ios_base::app);
  textFile << numAtoms << std::endl;

#endif

  for (int j=0;j<numAtoms;j++){
    //
    //    write out phase space coordinates
    //
    data_file << posTrial[j].x << posTrial[j].y << posTrial[j].z;
    data_file << vel[j].x << vel[j].y << vel[j].z;
#ifdef TEXTFILE
    textFile << posTrial[j].x << " " << posTrial[j].y << " " << posTrial[j].z << std::endl;
#endif
  }

  data_file.flush();
  data_file.close();
  textFile.close();
  writeConfigurationVMD(index);
  wroteImage[index] = true;


}


void ProteinChain::runMC(int steps){
  int accepts = 0;
  

  for (int i=0;i<numAtoms;i++) posTrial[i] = pos[i];
  potEnergy = computeTrialPotential();
  //std::cout << "   Starting " << steps << " MC moves with pot = " << potEnergy << std::endl;
 
  for (int i=0;i<steps;i++){
    accepts +=  crankshaftMC();
    writeFrame();
    writeMovieVMD();
  }
  //std::cout << "  MC acceptance ratio = " << double(accepts)/double(steps) << std::endl;
  MDinitialized = false;
}

int ProteinChain::crankshaftMC(){

  long int configIndex = configuration(pos);
  double Sb = 0.0;
  if (biasFlag) Sb = Sbias[configIndex];

  for (int i=0;i<numAtoms;i++) posTrial[i] = pos[i];
  int index = 1 + int( (numAtoms-1) * rnd() );
  Vector base = pos[index];
  Vector z = pos[index]-pos[index-1];
  Vector zhat = z/z.nrm();

  double phi = 2.0*M_PI*rnd();
  Vector axis = phi*zhat;
  Matrix rot = Rodrigues(axis);
  for (int j=index+1; j < numAtoms; j++){
      
    Vector dr_old = pos[j]-base;
    Vector dr_new = rot*dr_old;
    posTrial[j] = base+dr_new;
  }
  long int trialConfigIndex = -1;
  double potE = computeTrialPotential(posTrial, trialConfigIndex);
 
  if (potE == std::numeric_limits<double>::max() ) return 0;
  if ( stateActive[trialConfigIndex] == false) return 0;
  double bde = potE - Sbias[currentState];
  
  double Strial = Sbias[trialConfigIndex];
  bde = Strial - Sbias[currentState];


  


  int retval = 0;
  double T;
   if ( bde > 14.0) {
     T=0.0;
   } else if (bde < -20) {
     T=1.0; //accept
   } else {
     T=exp(-bde); 
   }

  if ( rnd() <= T ) {
      for (int i=0; i < numAtoms; i++) pos[i] = posTrial[i];
      
      if (currentState != trialConfigIndex) numBondMC++;
      if (levelTarget[currentState] == false) previousState = currentState;
      

  #ifdef PRINT_MC
      std::cout << "  ** MC move accepted: " << currentState << " -> " << trialConfigIndex
          << " bde = " << bde << " Sbias is now " << Sbias[trialConfigIndex] << " index = " 
            << index << " angle = " << phi << std::endl;
  #endif
      currentState = trialConfigIndex;
      retval++;
    }
    else
    {
 #ifdef PRINT_MC
      std::cout << "  ----- MC move rejected for attempted move " << currentState << " -> " << trialConfigIndex 
          << " bde = " << bde << " Sbias is still " << Sbias[currentState]<< " atom index = " 
            << index << " angle = " << phi << std::endl;
  #endif
    }

   //std::cout << "   pE = " << potEnergy << "   retval = " << retval << " numBondMC = " << numBondMC << std::endl;
   bondsInitialized = false; // don't use dynamical info
   return retval;
}

void ProteinChain::testCrankshaft(){

  int index = int(numAtoms/2);
  for (int i=0;i<numAtoms;i++){
    posTrial[i] = pos[i];
  }

  Vector base = pos[index];
  Vector z = pos[index]-pos[index-1];
  Vector zhat = z/z.nrm();
  
  int nrot = 100;
  double dphi = 2.0*M_PI/double(nrot);
  for (int i=1; i<= nrot; i++){
    double phi = dphi*double(i);
    Vector axis = phi*zhat;
    Matrix rot = Rodrigues(axis);

    // double det = rot.det();
    //std::cout << "  Determinant of rotation matrix is " << det << std::endl;
    for (int j=index+1; j < numAtoms; j++){
      
      Vector dr_old = posTrial[j]-base;
      Vector dr_new = rot*dr_old;
      pos[j] = base+dr_new;
      //std::cout << "  position " << j << " is now " << pos[j] << std::endl;
    }
    writeFrame();
    writeMovieVMD();
    
  }

  for (int i=0; i < numAtoms; i++) pos[i] = posTrial[i];
  
 

}

void ProteinChain::resetSystem(){
  std::cout << "   Recovering after error at time " << timeNow << std::endl;
  timeNow = checkpointTime - checkpointInterval; // time of last checkpoint
  std::cout << "  Setting time back to " << timeNow << std::endl;

  std::exit(EXIT_FAILURE);
  loadConfiguration();
  initializeVelocities(); // use new velocities
  for (int i=0;i<numAtoms;i++) localTime[i] = timeNow;
  setupEvents();
  
  potEnergy = computePotential(pos, currentState);
  totEnergy1  = computeEnergy(pos, kineticEnergy1, potEnergy1);
  totEnergy = totEnergy1;
  initialEnergy = totEnergy1;
  kineticEnergy = kineticEnergy1;
  potEnergy = potEnergy1;

  std::cout << timeNow << " " << totEnergy1 << " " << kineticEnergy1 << " " << potEnergy1 << " " << potEnergy << std::endl;
}

bool ProteinChain::computeLevelEntropy(int& updateInterval){

    bool done = false;
    int ns = 2; // only 2 states
    int totalCounts = configurationCount[startingState]+configurationCount[endingState];
    std::cout << "  Qsize " << hybridQ->getAverageQSize() << " Level entropy counts for states " << startingState 
              << " - " << endingState << std::endl;
    double norm = 0.5*( configurationCount[startingState]+configurationCount[endingState] ); // expected count for uniform distribution
    double av1 = configurationCount[startingState]/norm;
    double av2 = configurationCount[endingState]/norm;
 
    std::cout << " starting state " << startingState << "   has count " << configurationCount[startingState] << "  ratio = " 
                  << configurationCount[startingState]/norm << " norm = " << norm << std::endl;
    std::cout << " ending state " << endingState << "   has count " << configurationCount[endingState] << "  ratio = " 
                  << configurationCount[endingState]/norm << " norm = " << norm << std::endl;
  
    //   Use G-test to test for uniformity
    double gval = 100000.0;
    if (configurationCount[startingState] > 0) gval  = configurationCount[startingState]*log(configurationCount[startingState]/norm);
    if (gval < 100000.0 and configurationCount[endingState] > 0) gval += configurationCount[endingState]*log(configurationCount[endingState]/norm);

    double q2 = 1.0 + (ns+ 1.)/(6.0*totalCounts) + double(ns*ns)/(6.0*totalCounts*totalCounts);
    //std::cout << " finite size correction q2 = " << q2 << std::endl; 
    gval *= 2.0/q2;

    double cdf = gsl_cdf_chisq_P(gval, 1);
    double pval = 1.0 - cdf;
    int sIndex = stateIndexToState[startingState];
    int nVisits = states[sIndex].getVisits(activeBond);
    


    std::cout << std::endl << "  Level entropy Computed pval = " << pval << " for cycleNumber " << cycleNumber << "  gval = " << gval
              << "  nvisits = " << nVisits << " for active bond " << activeBond 
              << " av1 = " << av1 << "  av2 = " << av2 << " scaler = " << entropyScaler << std::endl;

   //computeLevelPassageTimes();

    int numDistances = distances.size();
    bool distancesSatisfied = false;
    if (numDistances >= maxNumDistances)
    {
		  distancesSatisfied = true;
    }
    else
    {
      std::cout << "  Number of distances " << numDistances << " less than " << maxNumDistances << " so continue iterations." << std::endl;
    }

    //if (pval > pcrit and (numDistances > 0.75*maxNumDistances)  and (nVisits > configsPerState/3) ){
    
  if (distancesSatisfied and pval > pcrit )
  {
    if (iterateGamma and entropyScaler > minEntropyScaler)
    {
      done = false;
      entropyScaler *= 0.5; // reduce entropy scaling
      updateInterval *= 2;
      if (updateInterval > 800) updateInterval = 800;
    }
    else 
    {
      done = true;
      std::cout << "  Convergence found and entropy stored in entropy.dat.  Number of outer collisions is " << nVisits << std::endl;
      std::cout << "  Entropy difference: S" << startingState << " - S" << endingState << " = " << Sbias[startingState] - Sbias[endingState]  << std::endl;
      checkConvergenceFlag = false;
      doneScanIteration = true;
      if (verifyConvergence) {
          doneScanIteration = false;
          checkConvergenceFlag = true;
          done = false;
          std::cout << "  Will now verify the convergence with additional test." << std::endl;
      } else {

  #ifndef TIMING   
            computeLevelPassageTimes();
  #endif
            doneScanIteration = true; // finished with this iteration of starting state
            return done;
        }
        cycleNumber = 0;
      }

  } else {
    done = false;
    checkConvergenceFlag = false;
    std::cout << " Not converged yet." << std::endl;
  }


  // update Sbias
  norm = 0.0;
  for (int i=0;i<numStates;i++){
      norm += configurationCount[i];
  }
  norm /= 2.0;
  if (norm == 0.0){
    std::cerr << " No points! Setting norm to unity." << std::endl;
    norm = 1.0;
  }
  //
  //   Leave Sbias of starting state the same
  //
  if (configurationCount[endingState] > 0) {
    if (configurationCount[startingState]>0){
      Sbias[endingState] += entropyScaler*log( double(configurationCount[endingState])/double(configurationCount[startingState]));
    } else {
      Sbias[endingState] += entropyScaler*log(configurationCount[endingState]);
    }
  } else {
    Sbias[endingState] -= entropyScaler*log(norm);
  }

  // if (pval > 0.0001) 
  // { 
  //   entropyScaler *= 0.9;
  // } 
  // else 
  // {
  //   entropyScaler *= 1.05;
  // }
  // if (entropyScaler < minEntropyScaler) entropyScaler = minEntropyScaler;
  // if (entropyScaler > maxEntropyScaler) entropyScaler = maxEntropyScaler;

  //if (configurationCount[startingState]>0){
  //  Sbias[startingState] += entropyScaler*log(configurationCount[startingState]/norm);
  //} else {
  //  Sbias[endingState] -= entropyScaler*log(norm);
  //}
  std::cout << " Sbias[" << startingState << "] = " << Sbias[startingState] << std::endl;
  std::cout << " Sbias[" << endingState << "] = " << Sbias[endingState] 
        << "  delta S = " << Sbias[startingState]  - Sbias[endingState] << std::endl;

  std::ofstream entropyFile("entropy.dat");
  for (int i=0;i<numStates;i++){
      //
      std::bitset<numberBits> bondState(i);
      std::vector<bool> state;
      for (int j=0;j<numBonds;j++) state.push_back( (bool)bondState[j] );
      std::string bondString = bondState.to_string();
      double energy_i = -1.0*bondState.count();

      entropyFile << i << " " << energy_i << " " << setprecision(10) << Sbias[i]-Sbias[numStates-1] << " " << state << std::endl;
  }

  entropyFile.close();

  for (int i=0;i<numStates;i++)configurationCount[i] = 0;
  return done;
}


void ProteinChain::computeEntropy(){
  //
  // output entropy
  //
  double norm = 0.0;
  
  for (int i=0;i<numStates;i++) {
      if (allowedState[i]) norm += configurationCount[i];
  }
  int totalCounts = norm;
  norm /= double(numAllowedStates); // average occupation number for uniform sampling
  

  double S0 = 0.0;
  
  int firstZeroState = numStates;
  double eFirst = -numBonds;
  double gval = 0.0;
  double min = std::numeric_limits<double>::max();
  double Xisq = 0.0;
  //   Use G-test to test for uniformity
 
  int minState = 0;
 
  for (int i=0;i<numStates;i++){
    entropyMatrix[i][cycleNumber] = configurationCount[i];

    if ( configurationCount[i] < min and allowedState[i] ) {
      min = configurationCount[i];
      minState = i;	
    }
    
    if (configurationCount[i] > 0) {
        double gcontribution = 2.0*configurationCount[i]*log(configurationCount[i]/norm);
        gval += gcontribution;
#ifdef PRINT_GVALUE
        std::cout << "  State " << i << "  with count " << configurationCount[i] << " contributes " << gcontribution << " to g-value = " << gval << std::endl;
#endif
        double delta = configurationCount[i] - norm;	
        Xisq += delta*delta/norm;
    }
  } 
  double q2 = 1.0 + (numAllowedStates + 1.)/(6.0*totalCounts) + 
              double(numAllowedStates * numAllowedStates)/(6.0*totalCounts*totalCounts);
  //std::cout << " finite size correction q2 = " << q2 << std::endl; 
  gval /= q2; // sample size correction
  
  double pval = 1.0 - gsl_cdf_chisq_P(gval, numAllowedStates-1);
  double gcrit = gsl_cdf_chisq_Qinv(0.05, numAllowedStates-1);
  
  if (pval > 0.0000001) {
     entropyScaler *= 0.85;
    if (entropyScaler < 0.075) entropyScaler = 0.075;
  } else if (pval < 0.00000001) {
    entropyScaler *= 1.1;
    if (entropyScaler > 0.25) entropyScaler = 0.25;
    
  }
 
  std::cout << std::endl << std::endl << "  Computed pval = " << pval << " for cycleNumber " << cycleNumber << "  gval = " << gval
	      << " and Xi_sq = " << Xisq << " min fraction = " << min/norm << " for state " << minState << std::endl
	      << "  average number states visited = " << norm << " scaler = " << entropyScaler
	            << std::endl;
 
  cycleNumber++;
  if (cycleNumber >= cyclesPerAverage){
    cycleNumber = 0; 
    // update Sbias 
    vector<double> averageCount(numStates);

    norm = 0.0; 
    for (int i=0;i<numStates;i++){
      averageCount[i] = 0;
      for (int j=0;j<cyclesPerAverage;j++) averageCount[i] += entropyMatrix[i][j];
      norm += averageCount[i];
    }
    totalCounts = norm;
    //
    //   Expected number of counts for state i is Ei = (total number states sampled) * Pi = (total number states)/numStates = norm/numStates
    //   Observed number of counts for state i is Oi = averageCount[i]
    // 
    gval = 0.0;
    for (int i=0;i<numStates;i++) {
        if (!allowedState[i]) continue; // skip forbidden state

        double Ei = norm/numAllowedStates;
        double Oi = averageCount[i];
        if (Oi < 1) Oi = 1;
        double gcontrib = 2.0*Oi*log(Oi/Ei);
        gval += gcontrib;
#ifdef PRINT_GVALUE
        std::cout << "  State " << i << "  with count " << Oi << " contributes " << gcontrib << " to Global g-value " << gval << std::endl;
#endif
    }
    q2 = 1.0 + (numAllowedStates + 1.)/(6.0*totalCounts) + 
              double(numAllowedStates * numAllowedStates)/(6.0*totalCounts*totalCounts);
    //std::cout << " finite size correction q2 = " << q2 << std::endl; 
    gval /= q2;

    double pval_global = 1.0 - gsl_cdf_chisq_P(gval, numAllowedStates-1); 
    std::cout << "  Global pval = " << pval_global << " for gval = " << gval << std::endl;

    if (pval > pcrit or pval_global > pcrit) {
      std::cout << "  Convergence found and entropy stored in entropy.dat " << std::endl;

      if (numStates != numAllowedStates){
        reduceEntropy();
      }
      exit(0);
    }

    
    
    if (averageCount[zeroState] > 0){
      double O0 = averageCount[zeroState];
      double E0 = norm/numAllowedStates;
      S0 = Sbias[zeroState] + entropyScaler*log(O0/E0);
    } else {
      double E0 = norm/numAllowedStates;
      double dS = -log(E0);
      if (dS < -2.0) dS = -2.0; 
      S0 = Sbias[zeroState]+dS; 
    }

    
    std::ofstream entropyFile("entropy.dat");

    for (int i=0;i<numStates;i++){
  
      ///   	   Update bias
      double dS = 0.0;
      double Si_est = 0.0;
      if (averageCount[i] > 0) {
        double Oi = averageCount[i];
        double Ei = norm/numAllowedStates; 
        dS = entropyScaler*log(Oi/Ei);
	      Si_est = Sbias[i] + dS;
       
      } else {
        double E0 = norm/numAllowedStates;
        dS = -log(E0);
        if (dS < -2.0) dS = -2.0; 
      	Si_est = Sbias[i] + dS; // lower the entropy to make unobserved states more likely
	
      }	

      if (allowedState[i]) Sbias[i] = Si_est - S0;

      if (allowedState[i]) std::cout << "   Bias " << i << " is now " <<  Sbias[i] << "   Delta S = " << dS << "  ratio = " << averageCount[i]/norm*numAllowedStates << std::endl;
      //
      
      
      std::bitset<numberBits> bondState(i);
      std::vector<bool> state;
      for (int j=0;j<numBonds;j++) state.push_back( (bool)bondState[j] );
      std::string bondString = bondState.to_string();
      double energy_i = -1.0*bondState.count();

      entropyFile << i << " " << energy_i << " " << setprecision(10) << Sbias[i] << " " << i << " " << state << std::endl;
      
    }
  
  }

  for (int i=0;i<numStates;i++)configurationCount[i] = 0;
}

void ProteinChain::whamEntropy(){

  int maxIteration = 1000;
  double EPS = 1.0E-12;
 

  double Zbeta = 1.0;
  double *lambda = new double[numStates+1];
  double *lambdaOld = new double[numStates+1];
  double *energy = new double[numStates+1];


  int numPoints = 0;
  for (int i=0;i<numStates;i++) {
    lambdaOld[i] = 1.0;
    bitset<numberBits> integerIndex(i);
    int nb = integerIndex.count();
    energy[i] = -1.0*nb;
    numPoints += configurationCount[i];

    //std::cerr << " state " << i << " has count " << configurationCount[i]
    //<< std::endl;
  }

  double conv = 10000000.0;
  int iteration = 0;

  while ( (conv > EPS) and (iteration < maxIteration) ){

  
    Zbeta = 0.0;
    for (int k=0;k<numStates;k++) Zbeta += lambdaOld[k];
  
    for (int j=0;j<numStates;j++){
      lambda[j] = configurationCount[j]/double(numPoints)  * Zbeta;
    }
    conv = 0.0;
    double normalizer = lambda[numStates-1];

    for (int j=0;j<numStates;j++){
      if (normalizer > 0.0) lambda[j] /= normalizer;
      if (lambdaOld[j] > 0.0) conv += (1.0 - lambda[j]/lambdaOld[j])*(1.0 - lambda[j]/lambdaOld[j]);
    }

    iteration++;

    std::cout << " At iteration " << iteration << "  conv = " << conv << std::endl;

    for (int j=0;j<numStates;j++) lambdaOld[j] = lambda[j];


  }
  //
  // output entropy
  //
  double norm = double(numPoints);
  std::ofstream entropyFile("entropy.dat");
  if (configurationCount[numStates-1] > 0) norm = configurationCount[numStates-1];

  for (int i=0;i<numStates;i++){
    double Si = -1.0;
    double Si_est = -1.0;
    if (lambda[i] > 0.0) Si = log(lambda[i]) + Sbias[i];
    
    Si_est = Sbias[i] + log(configurationCount[i]/norm);

    std::bitset<numberBits> bondState(i);
    std::string bondString = bondState.to_string();

    double freq_i = configurationCount[i]/double(numPoints);

    entropyFile << i << " " << energy[i] << " " << Si-Sbias[numStates-1] << " " << bondString << " " 
		  << Si_est << " " << freq_i << "  " << configurationCount[i]/norm << std::endl;
  }


  delete[] lambda;
  delete[] lambdaOld;
  delete[] energy;


}


void remove(std::vector<double> &v)
{
        auto end = v.end();
        for (auto it = v.begin(); it != end; ++it) {
                end = std::remove(it + 1, end, *it);
        }

        v.erase(end, v.end());
}




bool ProteinChain::checkConvergence(){


   //   Use G-test to test for uniformity
    bool converged = false;
    int ns = 2; // only two states
    int totalCounts = configurationCount[startingState]+configurationCount[endingState];
    double norm = 0.5*( totalCounts); // expected count for uniform distribution
    for (int i=0;i<numStates;i++)  entropyMatrix[i][cycleNumber] = configurationCount[i];

    double gval = configurationCount[startingState]*log(configurationCount[startingState]/norm);
    gval += configurationCount[endingState]*log(configurationCount[endingState]/norm);
    double q2 =  1.0 + (ns + 1.)/(6.0*totalCounts) + ns*ns/(6.0*totalCounts*totalCounts); // finite size correction
    //std::cout << " finite size correction q2 = " << q2 << std::endl; 
    gval *= 2.0/q2;

    double cdf = gsl_cdf_chisq_P(gval, 1);
    double pval = 1.0 - cdf;

    gvalArray.push_back(gval);

    std::cout << "   ------  Cycle " << cycleNumber << " of " << cyclesPerAverage << " has g = " << gval << "  pval = " << pval << std::endl;
    std::cout << "       **** counts are " << configurationCount[startingState] << "   " << configurationCount[endingState] << std::endl;
    cycleNumber++;
    //computeLevelPassageTimes();
    
    if (cycleNumber >= cyclesPerAverage){
         int devIndex;

         std::sort(gvalArray.begin(), gvalArray.end());

         //  construct empirical cdf

         vector<double> ecdfArray;
         for (int i=0;i<gvalArray.size();i++){
            double value = gvalArray[i];
            int count = 0;
            for (int j=0;j<gvalArray.size();j++){
                if (gvalArray[j] < value) count++;
            }
            double ecdf_i = double(count)/gvalArray.size(); 
            ecdfArray.push_back(ecdf_i);
         }


         //remove(gvalArray); // remove duplicates
         
         int nv = gvalArray.size();
         for (int i=0;i<nv;i++){
             double cdf_i = gsl_cdf_chisq_P(gvalArray[i], 1);
             cdfArray.push_back(cdf_i);
         }

         double qval = pa_q(ecdfArray, cdfArray, devIndex);

         std::cout << "   q-value of gvalues is " << qval << " with max. deviation at " << devIndex << "  gval = " << gvalArray[devIndex] << std::endl;
         cycleNumber = 0;

        double *averageCount = new double[numStates];
        for (int i=0;i<numStates;i++) {
            averageCount[i] = 0.0;
            for (int j=0;j<cyclesPerAverage;j++) averageCount[i] += entropyMatrix[i][j];
        }


        std::ofstream entropyFile("entropy.dat");
        totalCounts = averageCount[startingState]+averageCount[endingState];

        norm = 0.5*( averageCount[startingState]+averageCount[endingState] ); // expected count for uniform distribution

        gval = averageCount[startingState]*log(averageCount[startingState]/norm);
        gval += averageCount[endingState]*log(averageCount[endingState]/norm);
        int ns = 2; // only 2 states
        double q2 = 1.0 + (ns + 1.)/(6.0*totalCounts) + double(ns*ns)/(6.0*totalCounts*totalCounts);
        gval *= 2.0/q2;

        cdf = gsl_cdf_chisq_P(gval, 1);
        double pval_global = 1.0 - cdf;
        std::cout << "   Global p-value = " << pval_global << std::endl;

        //double S0 = Sbias[startingState];
        for (int i=0;i<numStates;i++){

          ///          Update bias
             if (averageCount[i] > 0) {
                Sbias[i] = Sbias[i] + 0.2 * log( double(averageCount[i]) / double(averageCount[startingState]) );
             }
             if (averageCount[i] > 0) std::cout << "   Bias " << i << " is now " <<  Sbias[i] << "  ratio = " << averageCount[i] << std::endl;
             //


             std::bitset<numberBits> bondState(i);
             std::vector<bool> state;
             for (int j=0;j<numBonds;j++) state.push_back( (bool)bondState[j] );
             std::string bondString = bondState.to_string();
             double energy_i = -1.0*bondState.count();

             entropyFile << i << " " << energy_i << " " << setprecision(10) << Sbias[i]- Sbias[numStates-1] << " " << state << std::endl;

       }


        if (qval > 0.1 or pval_global > pcrit){
             std::cout << " Good fit." << std::endl;
             //for (int i=0;i<numStates;i++) std::cout << "  Entropy " << i << " = " << Sbias[i] << std::endl;
             double deltaS = Sbias[endingState] - Sbias[startingState];
             std::cout << " Delta S for states S[" << startingState << "] - S[" << endingState << "] = " << deltaS << std::endl;
             computeLevelPassageTimes();
             doneScanIteration = true; // finished with this iteration of starting state
         }

       gvalArray.clear();
       cdfArray.clear();

       delete[] averageCount;
     }

     for (int i=0;i<numStates;i++)configurationCount[i] = 0;

     return converged;
}




double ProteinChain::pa_q(vector<double> &ecdf, vector<double>  &cdf, int &j)
{
  double a, b, d, e, c, g, h, q;
  int i;
  int n = ecdf.size();

  b = sqrt(n);
  e = g = h = q = 0.0;
  c = 2;
  j = 0;
  for (i = 0; i < n; ++i) {
    d = fabs(ecdf[i]-e);
    //d = fabs( double(i)/n - e);
    std::cout << "  i = " << i << "  ecdf = " << ecdf[i] << "   cdf = " << e << std::endl;
    if (d > h) {
      j = i;
      h = d;
    }
    e = cdf[i];
    d= fabs(e - ecdf[i]);
    //d = fabs(e-double(i)/n);
    if (d > h) {
      j = i;
      h = d;
    }
  }
  d = b*h+0.12*h+0.11*h/b;
  a = -2*d*d;
  for(i = 1;i<= 200;++i){
    d = c*exp(a*i*i);
    q += d;
    if (fabs(d) <= g)
      break;
    c = -c;
    g = fabs(d)/1000;
  }
  if (i > 200) std::cout << "#pa_q: no convergence." << std::endl;
  //std::cerr << " in new K-S, h = " << h << " and q = " << q << std::endl;
  return q;
}


