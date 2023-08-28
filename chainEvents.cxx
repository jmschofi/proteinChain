#include "chainMol.h" 
//double debugTime = 3069.0;
//double debugTime = 0.0;
int debugId1 = 16;
int debugId2 = 17;

//#define PRINT_TRANSITIONS
//#define PRINT_EVENT
//#define PRINT_BONDING
//#define PRINT_ENTRIES
//#define PRINT_BOND_CHANGE
//#define PRINT_BONDING_CALC
//#define PRINT_SW
//#define PRINT
//#define TEST_NEW
//#define USE_OLDQ

void ProteinChain::setupEvents(){
    Vector R(0,0,0);
	for (int i=0;i<numAtoms;i++){
		double dt = timeNow - localTime[i];
		pos[i] += vel[i]*dt;
        R += pos[i];
		localTime[i] = timeNow;
		interactionNumber[i] = 0; // NEW
	}
    R /= double(numAtoms);

        // zero the center of mass
    for (int i=0;i<numAtoms;i++) pos[i] -= R;
#ifdef USE_OLDQ      
	InitEventList(numAtoms); // resetting system after position update - not really needed
#else
	hybridQ->reset();
#endif

	for (int na = 0; na < numAtoms ; na++ ) calculateAllCollisions(na, -2);
	calculateAllBondCollisions();

#ifdef USE_OLDQ
	ScheduleEvent(-1, -1, runTime, 0.0,  STOP, -1); 
#else
	hybridQ->ScheduleEvent(-1, -1, runTime, 0.0,  STOP, -1, -1, -1);
#endif
	
	if (configFlag) {
#ifdef USE_OLDQ
		ScheduleEvent( -1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1);// measurements and output have idA = -1
#else
		hybridQ->ScheduleEvent(-1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1, -1, -1);
#endif

		//std::cout << "   Configuration to be recorded at time " << timeNow+configInterval << " runtime = "<< runTime << std::endl;
	}
	if (dataFlag and distances.size() < maxNumDistances) {
#ifdef USE_OLDQ
		ScheduleEvent( -1, -1, timeNow + dataInterval, 0.0, DISTANCEMEASUREMENT, -1);// measurements and output have idA = -1
#else
		hybridQ->ScheduleEvent(-1, -1, timeNow + dataInterval, 0.0, DISTANCEMEASUREMENT, -1, -1, -1);
#endif
		

	}

#ifdef USE_OLDQ
	ScheduleEvent( -1, -1, timeNow + syncInterval, 0.0, SYNCHRONIZE, -1);
#else
	hybridQ->ScheduleEvent( -1, -1, timeNow + syncInterval, 0.0, SYNCHRONIZE, -1, -1, -1);
#endif
	//std::cout << " Set initial sync event for " << timeNow + syncInterval << std::endl;
	
	double pe, ke = 0.0;
	initialEnergy = computeEnergy(pos, ke, pe);
	MDinitialized = true;
	collider1 = collider2 = -1;

	//hybridQ->printQueue();
}

bool ProteinChain::processEvent(){
	bool done = false;
#ifdef USE_OLDQ
	DrawNextEventQ(); // places event in evSchedule structure
#else
	hybridQ->drawQueue(evSchedule, interactionNumber);

	double evTime = evSchedule.time;
	timeNow = evTime;
#endif

#ifdef PRINT_EVENT
 	std::cout << std::endl << std::endl << "  Event drawn at time " << timeNow << " is:" << std::endl;
	evSchedule.print();
#endif

#ifdef TEST_NEW
    std::cout << std::endl << std::endl << "  Event drawn at time " << timeNow << " is:" << std::endl;
	evSchedule.print();

	std::cout  << "  Test event from new queue is:" << std::endl;
	testEv.print();
	std::cout << std::endl;

	if (evSchedule.idA != testEv.idA and evSchedule.idB != testEv.idB and evSchedule.time != testEv.time){
		std::cerr << " Error in agreement on event from different queues." << std::endl;
                std::cout << std::endl << std::endl << "  Event drawn at time " << timeNow << " is:" << std::endl;
                evSchedule.print();

                std::cout  << "  Test event from new queue is:" << std::endl;
                testEv.print();
                hybridQ->printTopQueue();

		exit(1);
	}
#endif

	switch(evSchedule.collType){
		case STOP: 
		{
			
			for (int i=0; i < numAtoms; i++){
				double dt = timeNow - localTime[i]; 
				pos[i] += vel[i]*dt;
				localTime[i] = timeNow;
			}
			totEnergy1 = computeEnergy(pos, kineticEnergy1, potEnergy1);
			collider1 = collider2 = -1;
			done = true;
			//std::cout << " Energy at time " << timeNow << " is " << totEnergy1 << std::endl;
			break;
		} // end STOP case 
		case DISTANCEMEASUREMENT: 
		{
			//if ( (currentState == levelTarget1) or (currentState == levelTarget2))
			if (levelTarget[currentState] == true)
			{
				for (int i=0;i<numBonds;i++){
	                int i1 = bondIndex1[ i ];
        	        int i2 = bondIndex2[ i ];
               		Vector ri1 = pos[i1] + (timeNow - localTime[i1])*vel[i1];
                    Vector ri2 = pos[i2] + (timeNow - localTime[i2])*vel[i2];
                    Vector r12 = ri2 - ri1;
                    double r = r12.nrm();
					distanceMatrix[i].push_back(r);

				}

				int i1 = bondIndex1[ activeBond ];
				int i2 = bondIndex2[ activeBond ];
				Vector ri1 = pos[i1] + (timeNow - localTime[i1])*vel[i1];
				Vector ri2 = pos[i2] + (timeNow - localTime[i2])*vel[i2];
				Vector r12 = ri2 - ri1;
				double r = r12.nrm();

				//std::cerr << "  Recording bond distance " << r << " at time " << timeNow 
				//		<< " currentState = " << currentState << " for active bond " 
				//		<< activeBond << "  num recorded = " << distanceMatrix[0].size() << std::endl;

				distances.push_back( (float) r);
			}
			if (distances.size() < maxNumDistances){
#ifdef USE_OLDQ
					ScheduleEvent( -1, -1, timeNow + dataInterval, 0.0, DISTANCEMEASUREMENT, -1);// measurements and output have idA = -1	
#else
					hybridQ->ScheduleEvent(-1, -1, timeNow + dataInterval, 0.0, DISTANCEMEASUREMENT, -1, -1, -1);
#endif
			}

			//std::cout << " Set next distance measurement for " << timeNow + dataInterval << " current state " << currentState
				//	<<	" flag = " << levelTarget[currentState] << std::endl;
			break;
		}
		case SYNCHRONIZE: 
		{


			for (int i=0; i < numAtoms; i++){
				double dt = timeNow - localTime[i];
				pos[i] += vel[i]*dt;
				localTime[i] = timeNow;
			}

			// if (recordAllDistances)
			// {
			// 	for (int i=0;i<numAtoms;i++)
			// 	{
			// 		for (int j=i+2;j<numAtoms;j++)
			// 		{
						
			// 			double r = (pos[i]-pos[j]).nrm();
			// 			if (r < reachedMin[i][j]) {
			// 				reachedMin[i][j] = r; // just check with single rc
			// 				reachedMin[j][i] = r;
							
			// 			}
			// 		}
			// 	}
			// }

			if (timeNow >= outputTime){
				outputTime += outputInterval;
				totEnergy1 = computeEnergy(pos, kineticEnergy1, potEnergy1);
				numAv++;
				Eav += totEnergy1;
				PEav += potEnergy1;
				KEav += kineticEnergy1;
				printOutput(timeNow);
			}


			if (movieFlag) {
				if (timeNow >= movieTime){
					movieTime += movieInterval;
					writeFrame();
					writeMovieVMD();
				}
			}

			if (timeNow >= checkpointTime) {
				checkpoint();
				checkpointTime += checkpointInterval;
			}

			if (timeNow >= redrawVelocityTime){
				redrawVelocityTime += redrawVelocityInterval;
				if (redrawVelocityFlag) initializeVelocities();
			}

			//std::cout << " Executed sync event at time " << timeNow << "  next outputTime is " << outputTime << std::endl;
			setupEvents();
			break;


		} // end SYNCHRONIZE case 

		case CONFIGMEASUREMENT: 
		{

#ifdef USE_OLDQ
			ScheduleEvent( -1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1);// measurements and output have idA = -1
#else
			if (timeNow + configInterval <= runTime) hybridQ->ScheduleEvent(-1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1,-1, -1);
#endif

			if (generateImagesFlag and !noAppend) saveStructure(currentState);

			if (generateImagesFlag and noAppend and !wroteImage[currentState] ) saveStructure(currentState);

			localHist[currentState]++;

			//std::cout << "  recorded state " << currentState << " at time " << timeNow << " next measurement at " << timeNow + configInterval
			//	<< "  interval = " << configInterval << std::endl;
					
			break;

		} // end CONFIGEASUREMENT case

		case BOND_COLLISION: 
		{
			collider1 = evSchedule.idA;
			collider2 = evSchedule.idB;
			interactionNumber[collider1]++;
			interactionNumber[collider2]++;

			double dt1 = timeNow - localTime[collider1];
			double dt2 = timeNow - localTime[collider2];

			pos[collider1] += vel[collider1]*dt1;
			pos[collider2] += vel[collider2]*dt2;

			localTime[collider1] = timeNow;
			localTime[collider2] = timeNow;

			bitset<numberBits> stateIndex(currentState);
			bitset<numberBits> toIndex(currentState);
			toIndex.flip( evSchedule.bondIndex );
			long int toState = toIndex.to_ulong();

			int sIndex = stateIndexToState[currentState]; // index of sorted state
			//int sIndex = currentState;
			states[sIndex].incrementVisits(evSchedule.bondIndex); // increment count of rc collisions for state and bond

			int startState = currentState;

			double delta_S = Sbias[toState] - Sbias[currentState];
			double dU = delta_S;  // entropy used as potential

			

			// if ( (bondActive[ evSchedule.bondIndex ] == false) or (stateActive[toState]==false) )
			// {
			// 	dU = 10000000000.0; // bounce back
			// }

			if (stateActive[toState] == false) dU = std::numeric_limits<double>::max();
 
			double preKe = 0.5*mass*(vel[collider1].nrm2() + vel[collider2].nrm2() );
	  
			Vector r12=pos[collider2] - pos[collider1];
			Vector r12hat = r12/r12.nrm();
			Vector v12=vel[collider2] - vel[collider1];
			double Vr = r12hat*v12;
	  
			double discr = Vr*Vr - 4.0*dU/mass;
			double alpha = 0.0;

#ifdef PRINT_BONDING
			std::cerr << std::endl << "  Processing bond collision between " << collider1 << "-" << collider2
							<< " at time " << timeNow << " distance " << r12.nrm() << " e = "
							<< potEnergy << " mass = " << mass << std::endl << "  current state is " << currentState << "  to state is "
							<< toState << "  bond involved is " << evSchedule.bondIndex
							<< " and bonding state is " << bondExists[evSchedule.bondIndex]
							<< "=" << stateIndex[ evSchedule.bondIndex ] << std::endl
							<< "    from bitset " << stateIndex << "  to bitset " << toIndex << "   S_from = " << Sbias[currentState]
							<< "   S_to = " << Sbias[toState] << "    dS = " << delta_S << std::endl;
#endif

			if (discr < 0){
				// not enough kinetic energy to escape, so bounce back
				alpha = Vr;
				dU = 0.0;

		  
#ifdef PRINT_BONDING
				std::cout << "   **  Bounce back since discr = " << discr 
								<< "  Vr = " << Vr << " in well " << evSchedule.bondIndex
								<< "  current state is " << currentState << std::endl;
#endif
		  
				if ( !bondExists[ evSchedule.bondIndex ] ){
					if (evSchedule.deltaU < 0.0){
#ifdef PRINT_BONDING
						//  This scenario should be ok if the bonding state has higher entropy
						std::cout << "      Bounce back called from outside since entropy with additional bond is higher  Vr = "
										<< Vr << " dS = " << delta_S << std::endl;
						evSchedule.print();
#endif
					} else {
						std::cerr << "Error: bond should exist but doesn't dynamically:  Vr = " << Vr << std::endl;
						evSchedule.print();
						errorFlag = true;
						exit(1);
#ifdef USE_OLDQ
						ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
						hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
					}
				}

			} else { // discr > 0

				if (Vr < 0.0){
					alpha = 0.5*( Vr + sqrt(discr) ); // entering well
			  
					if ( bondExists[ evSchedule.bondIndex ] ){
						std::cout << " Error: Entering well called but bond already exists according to dynamical accounting."
									<< std::endl;
						std::cout << "   Vr = " << Vr << "   discr = " << discr << "  d = " << r12.nrm()
									<< "  beta = " << beta << std::endl;
						evSchedule.print();
						errorFlag = true;
						exit(1);
#ifdef USE_OLDQ
						ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
						hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
					}// bondExists

#ifdef PRINT_ENTRIES
					std::cout << "*** Entry into well " << evSchedule.bondIndex << "  Vr = " << Vr 
							<< " sindex = " << sIndex << " currentState = " << currentState << " dU = " << dU 
							<< " endingState = " << endingState << " discr = " << discr << " numBondEvents = " << numBondEvents << std::endl;
#endif

#ifdef PRINT_BOND_CHANGE
					std::cout << "*** Entry into well " << evSchedule.bondIndex << "  Vr = " << Vr 
							<< " currentState = " << currentState 
							<< " endingState = " << endingState << " dU = " << dU << " discr = " << discr << std::endl;
#endif
			  
					bondExists[ evSchedule.bondIndex ] = true;
					if (levelTarget[currentState] == false) previousState = currentState; // record last non-target state
					//
					currentState = toState;
					numBondEvents++;

					// if (levelTarget[currentState] )
					// {
					// 	if ( currentState == endingState )
					// 	{
					// 		// reached last active state, so previous state can reach 
					// 		//if (stateConnected[previousState] == 0) std::cout << "  State " << previousState << " can reach ultimate target " << endingState << std::endl;
					// 		stateConnected[previousState]++; // previous state is connected

					// 		if (dihedralOutput and (stateConnected[previousState] < 50000) )
					// 		{
					// 			for (int i=0;i<numAtoms;i++)
					// 			{
					// 				for (int i=0; i < numAtoms; i++)
					// 				{
					// 					double dt = timeNow - localTime[i];
					// 					posTrial[i] = pos[i] + vel[i]*dt;	
					// 				}

					// 			}

					// 			DihedralState dihedrals(numAtoms, posTrial);

					// 			//std::cout << "  Writing out dihedrals for state " << currentState << " at time " << timeNow 
					// 			//			<< " from previous originiating state " << previousState << std::endl;

					// 			std::bitset<numberBits> startState(startingState);
					// 			std::bitset<numberBits> levelState(previousState);
					// 			std::bitset<numberBits> difference = (startState ^= levelState); // XOR bit
					// 			int sindex = -1;
					// 			for (int i=0;i<numberBits;i++) 
					// 			{
					// 				if (difference[i]) sindex = i; // bit set is difference in bonding pattern
					// 			}
					// 			dihedralFile << sindex << " " << previousState << " ";
					// 			dihedrals.output(dihedralFile);
					// 		}
							
					// 	}
					// }
	

              	} else { // Vr > 0.0

                  	alpha = 0.50*( Vr - sqrt(discr) );
		  
					if ( !bondExists[ evSchedule.bondIndex ] ){
						std::cout << " Error: Leave well called but bond shouldn't exist according to dynamical accounting."
								<< std::endl;
						evSchedule.print();
						errorFlag = true;
						exit(1);
#ifdef USE_OLDQ
						ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
						hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
					} // bondExists
			  
					bondExists[ evSchedule.bondIndex ] = false;
					if (levelTarget[currentState] == false) previousState = currentState; // record last non-target state
					currentState = toState;
					
					numBondEvents++;

#ifdef PRINT_BOND_CHANGE
					std::cout << " $$$$  Leaving well " << evSchedule.bondIndex << "  Vr = " << Vr << " currentState = "
								<< currentState << " dU = " << dU << " discr = " << discr << std::endl;
#endif

            
				} // Vr > 0.0

			} // discr > 0.0
 
			potEnergy = Sbias[currentState];

			Vector dv = r12hat*alpha;
			vel[collider1] += dv;
			vel[collider2] -= dv;    
#ifdef PRINT_TRANSITIONS
			if (startState != currentState)
			{
				std::cout << " Transition from state with previous (" << previousState << ")   transition " << startState 
						<< " -> " << currentState << "  t = " << timeNow << std::endl;

				//if (currentState == endingState) exit(1);
			}
#endif

#ifdef PRINT_BONDING
			double postKe = 0.5*mass*(vel[collider1].nrm2() + vel[collider2].nrm2() );
			double Vrnew =  (vel[collider2] - vel[collider1])*r12hat;
			std::cout << "   E before collision = " << preKe << "    after collision " << postKe + dU
						<< "  potE = " << potEnergy << "  current state is " << currentState << std::endl;
	  
			if (fabs(postKe+dU - preKe)/preKe > 0.000000000001) {
				std::cerr << "  Error in energy conservation?" << std::endl;
				std::cerr << "  dU set to " << dU << " and evSchedule is " << std::endl;
				evSchedule.print();
				exit(1);
			}
#endif

			calculateBondCollisions(collider1, collider2);
			calculateAllCollisions(collider1, collider2);
			collider1 = collider2 = -1;
			break;
		} // end bond collision case
		case MIN_COLLISION: 
		case MAX_COLLISION: { // GENERAL HARD CORE COLLISIONS

			collider1 = evSchedule.idA;
			collider2 = evSchedule.idB;
			interactionNumber[collider1]++;
			interactionNumber[collider2]++;
			double dt1 = timeNow - localTime[collider1];
			double dt2 = timeNow - localTime[collider2];
    
			pos[collider1] += vel[collider1]*dt1;
			pos[collider2] += vel[collider2]*dt2;

			localTime[collider1] = timeNow;
			localTime[collider2] = timeNow;

			Vector r12=pos[collider2] - pos[collider1];
			Vector r12hat = r12/r12.nrm();
			Vector v12=vel[collider2] - vel[collider1];
			double b = r12hat*v12;
			Vector dv = r12hat*b;
#ifdef PRINT

			std::cout << "  When processing non-bonded collision, dijsq = " 
							<< r12.nrm2() << "  dv = " << dv <<  " b = " << b << std::endl;
			std::cout << "           collider1 = " << collider1 << "   bondedAtom = " << isBondedAtom[collider1] << std::endl;
			std::cout << "           collider2 = " << collider2 << "   bondedAtom = " << isBondedAtom[collider2] << std::endl;
#endif

			vel[collider1] = vel[collider1] + dv;
			vel[collider2] = vel[collider2] - dv;

			calculateAllCollisions(collider1, collider2);
			if ( isBondedAtom[collider1] or isBondedAtom[collider2] )  calculateBondCollisions(collider1, collider2); // atoms in bonded interactions

			collider2 = -1;
			collider1 = -1;

		}
	} // end of switch statements

	return done;
}

bool ProteinChain::processEventOld(){


  bool done = false;
  DrawNextEventQ(); // places event in evSchedule structure

  double evTime = evSchedule.time;
  timeNow = evTime;

   Event testEv;
   hybridQ->drawQueue(testEv, interactionNumber);

#ifdef TEST_NEW
    std::cout << "  Event drawn at time " << timeNow << " is:" << std::endl;
	evSchedule.print();

	std::cout  << "  Test event from new queue is:" << std::endl;
	testEv.print();
	std::cout << std::endl;
#endif

  if (dataFlag){
    if (timeNow >= dataTime){
      std::cerr << "   Data recording time is " << dataTime << "   timeNow = " << timeNow << std::endl;
      dataTime += dataInterval;
      writeData();
    }
  }

  if (evSchedule.collType == STOP){
    //std::cerr << "Executing STOP event at time " << timeNow << std::endl;
    for (int i=0; i < numAtoms; i++){
      double dt = timeNow - localTime[i];
      pos[i] += vel[i]*dt;
      localTime[i] = timeNow;
    }
    totEnergy1 = computeEnergy(pos, kineticEnergy1, potEnergy1);
    collider1 = collider2 = -1;
    //printOutput();
    done = true;
  } else if (evSchedule.collType == OUTPUT){

    for (int i=0; i < numAtoms; i++){
      double dt = timeNow - localTime[i];
      pos[i] += vel[i]*dt;
      localTime[i] = timeNow;
    }
    totEnergy1 = computeEnergy(pos, kineticEnergy1, potEnergy1);
    numAv++;
    Eav += totEnergy1;
    PEav += potEnergy1;
    KEav += kineticEnergy1;
    printOutput(timeNow);

   
    if (timeNow >= checkpointTime) {
      checkpoint();
      checkpointTime += checkpointInterval;
    }

    if (timeNow >= acfCorrelationTime){
      if (acfFlag) outputAcf();
    }

    if (timeNow >= redrawVelocityTime){
      redrawVelocityTime += redrawVelocityInterval;
      if (redrawVelocityFlag) initializeVelocities();

    }

    if (movieFlag) {
        if (timeNow >= movieTime){
            movieTime += movieInterval;
	    writeFrame();
            writeMovieVMD();
        }
    }
    collider1 = collider2 = -1;
    InitEventList(numAtoms); // resetting system after position update - not really needed
    setupEvents();
#ifdef USE_OLDQ
    ScheduleEvent(-1, -1, runTime, 0.0,  STOP, -1);
#else
	hybridQ->ScheduleEvent(1, -1, runTime, 0.0,  STOP, -1,-1,-1);
#endif
   
  } else if (evSchedule.collType == BOND_COLLISION) {

	  
	  collider1 = evSchedule.idA;
	  collider2 = evSchedule.idB;

	  interactionNumber[collider1]++;
	  interactionNumber[collider2]++;

	  double dt1 = timeNow - localTime[collider1];
	  double dt2 = timeNow - localTime[collider2];

	  pos[collider1] += vel[collider1]*dt1;
	  pos[collider2] += vel[collider2]*dt2;

	  localTime[collider1] = timeNow;
	  localTime[collider2] = timeNow;

	  bitset<numberBits> stateIndex(currentState);
	  bitset<numberBits> toIndex(currentState);
	  toIndex.flip( evSchedule.bondIndex );
	  long int toState = toIndex.to_ulong();

	  double delta_S = Sbias[toState] - Sbias[currentState];

#ifdef PRINT_BONDING
		  std::cerr << std::endl << "  Processing bond collision between " << collider1 << "-" << collider2
					<< " at time " << timeNow << " e = "
					<< potEnergy << " mass = " << mass << std::endl << "  current state is " << currentState << "  to state is "
					<< toState << "  bond involved is " << evSchedule.bondIndex
					<< " and bonding state is " << bondExists[evSchedule.bondIndex]
					<< "=" << stateIndex[ evSchedule.bondIndex ] << std::endl
					<< "    from bitset " << stateIndex << "  to bitset " << toIndex << "   S_from = " << Sbias[currentState]
					<< "   S_to = " << Sbias[toState] << "    dS = " << delta_S << std::endl;
#endif
 
	  double preKe = 0.5*mass*(vel[collider1].nrm2() + vel[collider2].nrm2() );
	  
	  Vector r12=pos[collider2] - pos[collider1];
	  Vector r12hat = r12/r12.nrm();
	  Vector v12=vel[collider2] - vel[collider1];
	  double Vr = r12hat*v12;
	  
	  //double dU = -1.0;
	  double dU = delta_S;
	  //if (Vr > 0.0) dU = 1.0;
	  //double dU = evSchedule.deltaU;


	  double discr = Vr*Vr - 4.0*dU/mass;
	  double alpha = 0.0;
   
	  if (discr < 0){
		  // not enough kinetic energy to escape, so bounce back
		  alpha = Vr;
		  dU = 0.0;
		  
#ifdef PRINT_BONDING
			  std::cout << "   **  Bounce back since discr = " << discr 
						<< "  Vr = " << Vr << " in well " << evSchedule.bondIndex
						<< "  current state is " << currentState << std::endl;
#endif
		  
		  if ( !bondExists[ evSchedule.bondIndex ] ){
		    if (evSchedule.deltaU < 0.0){
#ifdef PRINT_BONDING
				 //  This scenario should be ok if the bonding state has higher entropy
				std::cout << "      Bounce back called from outside since entropy with additional bond is higher  Vr = "
							<< Vr << " dS = " << delta_S << std::endl;
				evSchedule.print();
#endif
			} else {
		        std::cerr << "Error: bond should exist but doesn't dynamically:  Vr = " << Vr << std::endl;
		        evSchedule.print();
	         	errorFlag = true;
				exit(1);
#ifdef USE_OLDQ
				ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
				hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1,-1,-1);
#endif
			}
		}

	  } else {
        if (Vr < 0.0){
            alpha = 0.5*( Vr + sqrt(discr) ); // entering well
			  
            if ( bondExists[ evSchedule.bondIndex ] ){
				std::cout << " Error: Entering well called but bond already exists according to dynamical accounting."
								<< std::endl;
		     	std::cout << "   Vr = " << Vr << "   discr = " << discr << "  d = " << r12.nrm()
								<< "  beta = " << beta << std::endl;
                evSchedule.print();
                errorFlag = true;
#ifdef USE_OLDQ
                ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
				hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1,-1,-1);
#endif
            }// bondExists
			  
		  bondExists[ evSchedule.bondIndex ] = true;
		  if (currentState == startingState) numBondEvents++; // record entry out of most unbonded
		  if (levelTarget[currentState] == false) previousState = currentState; // record last non-target state
		  currentState = toState;
		  

#ifdef PRINT_BONDING
		  std::cout << "*** Entry into well " << evSchedule.bondIndex << "  Vr = " << Vr << std::endl;
#endif
          
        } else { // Vr > 0.0
			alpha = 0.50*( Vr - sqrt(discr) );
		  	if ( !bondExists[ evSchedule.bondIndex ] ){
		      	std::cout << " Error: Leave well called but bond shouldn't exist according to dynamical accounting."
					<< std::endl;
	        	evSchedule.print();
		    	errorFlag = true;
#ifdef USE_OLDQ
		    	ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
				hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1,-1,-1);
#endif

	        } // bondExistst
			  
		  	bondExists[ evSchedule.bondIndex ] = false;
			if (levelTarget[currentState] == false) previousState = currentState; // record last non-target state
		  	currentState = toState;

#ifdef PRINT_BONDING
		    std::cout << " $$$$  Leaving well " << evSchedule.bondIndex << "  Vr = " << Vr << std::endl;
#endif
            
        } // Vr > 0.0
	    potEnergy += dU;
	} // discr > 0.0

	Vector dv = r12hat*alpha;

	vel[collider1] += dv;
	vel[collider2] -= dv;    

#ifdef PRINT_BONDING
    double postKe = 0.5*mass*(vel[collider1].nrm2() + vel[collider2].nrm2() );
	double Vrnew =  (vel[collider2] - vel[collider1])*r12hat;
	std::cout << "   E before collision = " << preKe << "    after collision " << postKe + dU
	                << "  potE = " << potEnergy << "  current state is " << currentState << std::endl;
	  
	if (fabs(postKe+dU - preKe) > 0.000000000001) {
		std::cerr << "  Error in energy conservation?" << std::endl;
		std::cerr << "  dU set to " << dU << " and evSchedule is " << std::endl;
		evSchedule.print();
		exit(1);
	}
#endif

	calculateBondCollisions(collider1, collider2);
	calculateAllCollisions(collider1, collider2);
	collider1 = collider2 = -1;

} else if (evSchedule.collType == CONFIGMEASUREMENT){

#ifdef USE_OLDQ
	 ScheduleEvent( -1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1);// measurements and output have idA = -1
#else
	  hybridQ->ScheduleEvent(-1, -1, timeNow + configInterval, 0.0, CONFIGMEASUREMENT, -1,-1, -1);
#endif
	  localHist[currentState]++;
	  //std::cerr << "  Recording configuration " << currentState << " at time " << timeNow << std::endl;

  } else if (evSchedule.collType == ACFMEASUREMENT){
    
	// std::cout << "!!!  measurement of acf at " << timeNow << std::endl;
#ifdef USE_OLDQ
	ScheduleEvent( -1, -1, timeNow + acfInterval, 0.0, ACFMEASUREMENT, -1);// measurements and output have idA = -1
#else
	hybridQ->ScheduleEvent(1, -1, timeNow + acfInterval, 0.0, ACFMEASUREMENT, -1, -1, -1);
#endif

    int indexDifference = 30;
    int index1 = int(numAtoms/2) - indexDifference;

    for (int i = index1; i < index1+indexDifference; i++){
      int j = i + indexDifference;
      Vector rij = pos[j] - pos[i] + (timeNow-localTime[j])*vel[j] - (timeNow - localTime[i])*vel[i]; // update positions
      double rij_nrm = rij.nrm();
      Vector rijhat = rij/rij_nrm;
      double vij_rel = (vel[j]-vel[i])*rijhat;
      vij_rel_av += vij_rel;
      acf->addpoint(vij_rel);
      acfPoints++;
      //std::cout << "   acf point " << acfPoints << " added " << vij_rel << std::endl;
    }
    collider1 = collider2 = -1;

  } else {
    collider1 = evSchedule.idA;
    collider2 = evSchedule.idB;
    interactionNumber[collider1]++;
    interactionNumber[collider2]++;
    double dt1 = timeNow - localTime[collider1];
    double dt2 = timeNow - localTime[collider2];
    
    pos[collider1] += vel[collider1]*dt1;
    pos[collider2] += vel[collider2]*dt2;

    localTime[collider1] = timeNow;
    localTime[collider2] = timeNow;

    Vector r12=pos[collider2] - pos[collider1];

    Vector r12hat = r12/r12.nrm();
    Vector v12=vel[collider2] - vel[collider1];
    double b = r12hat*v12;
    Vector dv = r12hat*b;


    vel[collider1] = vel[collider1] + dv;
    vel[collider2] = vel[collider2] - dv;

    calculateAllCollisions(collider1, collider2);
    if ( isBondedAtom[collider1] or isBondedAtom[collider2] )  calculateBondCollisions(collider1, collider2); // atoms in bonded interactions

    collider2 = -1;
    collider1 = -1;
  }
  return done;
}




void ProteinChain::calculateAllBondCollisions(){
  //
  //   Positions updated so localTime = timeNow for all particles
  //
  double collisionTime = INFINITY;
  for (int j=0; j< numBonds; j++){
      int i1 = bondIndex1[j];
      int i2 = bondIndex2[j];
      Vector rij = pos[i2] - pos[i1];
      Vector vij = vel[i2] - vel[i1];

      double rijsq = rij.nrm2();
      double vijsq = vij.nrm2();

      double alpha = (rij*vij) /vijsq; // dot product
      double term = (rijsq - bondDistance[j]*bondDistance[j])/vijsq;
      double discr = alpha*alpha - term;

#ifdef PRINT_BONDING_CALC
      std::cerr << "  Computation for Bond " << j << "  between " << i1 << "-" << i2 << " at distance " << bondDistance[j] << " has bond state " << bondExists[j]
	           << " rijsq = " << rijsq << "  alpha = " << alpha << " discr = " << discr 
		    << "  times " << localTime[i1] << "=" << localTime[i2] << "=" << timeNow << std::endl;

#endif

      if (discr > 0.0){
		  //  two roots exist, otherwise none
		  double part1 = -alpha-sqrt(discr);
		  double part2 = -alpha+sqrt(discr);

		  if ( ((i1 == collider1) and (i2==collider2)) or ( (i1==collider2) and (i2==collider1)) ){
			if (alpha >= 0.0 and (!bondExists[j]) and (evSchedule.bondIndex == j) ) continue; // no collision possible as moving out of well
				if ( evSchedule.bondIndex == j) {
                	collisionTime = part2;
            	} else {
                    collisionTime = part1;
                }
#ifdef PRINT_BONDING_CALC
		        std::cout << "  Bonding calc. for leaving event " << i1 << "-" << i2
				   << " bond number " << j << " rijsq = " << rijsq << " alpha = " << alpha 
				   << "  term = " << term
				   << "  discr = " << discr << "  part1 (should be almost zero) = "
		   	           << part1 << "  part2 = " << part2
			           << "  coll. time = " << collisionTime << std::endl;

#endif

#ifdef USE_OLDQ
			  ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j); // time pair will leave well
#else
			  hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
		  } else {
			  if ( (term > 0.0) and (alpha < 0.0) ){
				  collisionTime = part1; // outside of well, moving into well at time t=part1
#ifdef PRINT_BONDING_CALC
			    std::cout << "  Bonding calc. for entering event " << i1 << "-" << i2
				   << " alpha = " << alpha << "  term = " << term
				   << "  discr = " << discr << "  part1 (positive) = "
		     	           << part1 << "  part2 (positive) = " << part2
				   << "  coll. time = " << collisionTime << std::endl;

#endif
#ifdef USE_OLDQ
				  ScheduleEvent(i1, i2, timeNow+collisionTime, -1.0,  BOND_COLLISION, j);
#else
				  hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, -1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
			  } else if (term < 0.0) {
				  // already inside a bonding well, 2 roots: part1<0.0 (entry time in past) and part2>0 (well exit)
				  collisionTime = part2;
#ifdef PRINT_BONDING_CALC
				 std::cout << "  Bonding calc. for outer event from inside: "
					   << i1 << "-" << i2 << " alpha = " << alpha 
					   << "  term = " << term
					   << "  discr = " << discr << std::endl << "      part1 (negative) = "
					   << part1 << "  part2 (positive) = " << part2
					   << "  coll. time = " << collisionTime << std::endl;

#endif
#ifdef USE_OLDQ
				  ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j);
#else
				  hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
			  }
		  }
      }
  } // sum over bonds

}

void ProteinChain::calculateBondCollisions(int ind1, int ind2){
  //
  //    called after a collision, so not all particle positions are current
  //
  double collisionTime = INFINITY;
  for (int j=0; j < numBonds; j++){
	int i1 = bondIndex1[j];
	int i2 = bondIndex2[j];

	if ( (i1 == ind1) or (i2 == ind1) or (i1 == ind2) or (i2 == ind2) ) {
			
		double dt1 = timeNow - localTime[i1];
		double dt2 = timeNow - localTime[i2];
		Vector drij = vel[i2]*dt2 - vel[i1]*dt1;
      
		Vector rij = pos[i2] - pos[i1] + drij;
		Vector vij = vel[i2] - vel[i1];

		double rijsq = rij.nrm2();
		double vijsq = vij.nrm2();
          

		double alpha= (rij*vij) /vijsq; // dot product
		double term = (rijsq - bondDistance[j]*bondDistance[j])/vijsq;            
		double discr = alpha*alpha - term;
			
#ifdef PRINT_BONDING_CALC
		std::cout << "  bond calculation " << j << " for pair " 
			  << i1 << "-" << i2 << " from calculationBondCollisions(" 
			  << ind1 << "," << ind2 << ") current bond state = "
			  << bondExists[j] << " term = " << term << "  alpha = "
			  << alpha << "  discr = " << discr << std::endl;
		if (discr > 0){
			double root1 = -alpha - sqrt(discr);
			double root2 = -alpha + sqrt(discr);
			std::cout << "        two roots are: " << root1 << "  and " << root2 << std::endl;
		}
#endif

		//  Check dynamic bond state
		if ( (evSchedule.bondIndex == j) and ((i1 == collider1) and (i2==collider2)) or ( (i1==collider2) and (i2==collider1)) ){
			// current colliding pair, so should be right at discontinuity in potential
			double part1 = -alpha-sqrt(discr);
			double part2 = -alpha+sqrt(discr);
#ifdef PRINT_BONDING_CALC
			std::cout << "        roots of current colliding pair " << i1 << "-" << i2
				  <<  " in recalculation are " << part1 << " and " << part2 
				  << "   alpha = " << alpha << " term = " << term << "  rij = " << rij.nrm() << std::endl;
#endif
		} else {
			if (bondExists[j] and (rijsq > bondDistance[j]*bondDistance[j] ) ){
			  
				std::cerr << "  Expected a bond to exist but d_ij = " << rij.nrm() << std::endl;
				std::cout << "  bond calculation " 
					  << i1 << "-" << i2 << " from calculationBondCollisions(" 
					  << ind1 << "," << ind2 << ") current bond state = " << bondExists[j]
					  << "  collider1 = " << collider1 << "   collider2 = " << collider2
					  << " term = " << term << "  time = " << timeNow << std::endl;
				std::cout << "  Current state is " << currentState << std::endl;
				std:cout << " rijsq - bondDistance**2 = " << rijsq - bondDistance[j]*bondDistance[j]
					 << "  alpha (Vr) = " << alpha << std::endl;
				std::cout << "  roots are " << -alpha - sqrt(discr) << "  and " << -alpha + sqrt(discr) << std::endl;
				std::cout << "   Printing previous event:" << std::endl;
				evSchedule.print();
				if (errorCount > 0){
					errorCount = 0;
					errorFlag = true;
					std::cerr << "  Already had an error, so abandoning attempt." << std::endl;
#ifdef USE_OLDQ
					ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
					hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1,-1,-1);
#endif
					continue;
				} else {
					std::cerr << "   Trying to continue by scheduling an immediate bond collision..." << std::endl;
					errorCount++;
#ifdef USE_OLDQ
					ScheduleEvent(i1, i2, timeNow, 1.0,  BOND_COLLISION, j); 
#else
					hybridQ->ScheduleEvent(i1, i2, timeNow, 1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
					continue;
				}
			}
				
			if ( (!bondExists[j]) and (rijsq < bondDistance[j]*bondDistance[j]) ){
			    
				std::cerr <<   " Bond should exist but term = " << term << std::endl;
				std::cout << "  bond calculation " 
					  << i1 << "-" << i2 << " from calculationBondCollisions(" 
					  << ind1 << "," << ind2 << ") current bond state = " << bondExists[j]
					  << " term = " << term << " time = " << timeNow << std::endl;
                                std::cout << "  Current state is " << currentState << std::endl;
			        std::cout << " rijsq - bondDistance**2 = " << rijsq - bondDistance[j]*bondDistance[j]
				         << "  alpha (Vr) = " << alpha << "  rij = " << sqrt(rijsq) << " bond distance is "
						 << bondDistance[j] << std::endl;

				std::cout << "   Printing previous event:" << std::endl;
				evSchedule.print();
				exit(1);
				if (errorCount > 0){
					errorCount = 0;
					errorFlag = true;
					std::cerr << "  Already had an error, so abandoning attempt." << std::endl;
#ifdef USE_OLDQ
					ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
					hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
					continue;
				} else {
					std::cerr << "   Trying to continue by scheduling an immediate bond collision..." << std::endl;
					errorCount++;
#ifdef USE_OLDQ
					ScheduleEvent(i1, i2, timeNow, -1.0,  BOND_COLLISION, j); 
#else
					hybridQ->ScheduleEvent(i1, i2, timeNow, -1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
					continue;
				}

			}
		}

		if (discr > 0.0){
			//  two roots exist, otherwise none
			double part1 = -alpha-sqrt(discr);
			double part2 = -alpha+sqrt(discr);

			if ( ((i1 == collider1) and (i2==collider2)) or ( (i1==collider2) and (i2==collider1)) ){

				double dU = 0.0;

#ifdef PRINT_BONDING_CALC
				std::cout << "  Pair bond calc. event for previously colliding pair "
					 << i1 << "-" << i2 
				        << "  bond index is " << j << "  last bond event was "
						<< evSchedule.bondIndex << " alpha = " << alpha << "  term = " << term
				 	<< "  discr = " << discr << std::endl << "          part1 = " << part1
				 	<< "  part2 = " << part2
				 	<< "  coll. time = " << part2  << " time = " << timeNow << std::endl;
#endif


				if (alpha >= 0.0 and (!bondExists[j]) ) {
#ifdef PRINT_BONDING_CALC
					std::cout << "  Not scheduling since moving away from bonding distance for bond "
							  << i1 << "-" << i2
							  << " alpha = " << alpha << "  term = " << term
							  << "  discr = " << discr << std::endl << "          part1 = " << part1
							  << "  part2 = " << part2
							  << "  coll. time = " << part2  << " time = " << timeNow << std::endl;
                                                 
#endif
		 			continue; // no collision possible as moving out of well
				}
                if (j == evSchedule.bondIndex){
                    // part1 = 0 since leaving current bonding well,  bondIndex is bond last processed
				    collisionTime = part2;
                } else {
                    //collisionTime = part1; // case where two beads have more than one bonding distance
                    //  part1 should be negative
					if ( j != evSchedule.bondIndex and part1 > 0.0) {
						collisionTime = part1; // entry rather than exit
						dU = -1.0;
					} else {
						collisionTime = part2;
						dU = 1.0; // exit
					}
#ifdef PRINT_BONDING_CALC
					std::cout << "  Staircase: previous collider inner collision time " << part1 
						<< "  outer collision time " << part2 << "  collision at "  << collisionTime << std::endl;
#endif
                    	
                }

                if (collisionTime < 0.0){
                    std::cerr << " Error!  negative collision time for bond calculation " << j << " scheduled last collision "
                                << evSchedule.bondIndex << "  root 1 = " << part1 << "  root 2 = " << part2 << std::endl;
                    exit(1);
                }

				if (!bondExists[j] and (evSchedule.bondIndex == j) ){
					std::cerr << " Computed leaving time for bond " << j << " that doesn't exist dynamically!" << std::endl;
					std::cerr << "  Pair bond calc. leaving event " << i1 << "-" << i2 
						  << " alpha = " << alpha << "  term = " << term
						  << "  discr = " << discr << std::endl << "         part1 (should be almost zero) = "
						  << part1 << "  part2 = " << part2
						  << "  coll. time = " << collisionTime  << " time = " << timeNow << std::endl;
					
					errorFlag = true;
#ifdef USE_OLDQ
					ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
					hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
				}
#ifdef USE_OLDQ
				ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j); // time pair will leave well
#else
				hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, dU,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif			
			} else {
				if ( (term > 0.0) and (alpha < 0.0) ){
					collisionTime = part1; // outside of well, moving into well at time t=part1
#ifdef PRINT_BONDING_CALC
						std::cout << "  Pair bond calc. for enter event "
							<< i1 << "-" << i2 << " alpha = " << alpha << "  term = " << term
							<< "  discr = " << discr << std::endl << "         part1 (pos) = "
							<< part1 << "  part2 (pos) = " << part2
							<< "  coll. time = " << collisionTime
							<< " time = " << timeNow << std::endl;

#endif

					if (bondExists[j]){
						   
						std::cerr << " Computed entry time for bond " << j << " that exists dynamically!" << std::endl;
						std::cerr << "  Pair bond calc. for enter event " << i1 << "-" << i2
							  << " alpha = " << alpha << "  term = " << term
							  << "  discr = " << discr << "  part1 (pos) = " << part1 << "  part2 (pos) = " << part2
							  << "  coll. time = " << collisionTime  << " time = " << timeNow << std::endl;
				                std::cout << "   Printing previous event:" << std::endl;
						evSchedule.print();
						errorFlag = true;
#ifdef USE_OLDQ
						ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
						hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
					}
#ifdef USE_OLDQ
					ScheduleEvent(i1, i2, timeNow+collisionTime, -1.0,  BOND_COLLISION, j);
#else
					hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, -1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
				} else if (term < 0.0) {
					// already inside a bonding well, 2 roots: part1<0.0 (entry time in past) and part2>0 (well exit)
					collisionTime = part2;
#ifdef PRINT_BONDING_CALC
					    std::cout << "  Pair Bond calc. for outer event from inside: "
						  << i1 << "-" << i2 << " for bond index " << j << " alpha = " 
						  << alpha << "  term = " << term
						  << "  discr = " << discr << std::endl
						  << "         part1 (neg) = "
						  << part1 << "  part2 (pos) = " << part2
						  << "  coll. time = " << collisionTime
						  << " time = " << timeNow << std::endl;
#endif
                                           
					if (!bondExists[j]){
						std::cerr << " Computed leaving time for bond " << j << " that doesn't exist dynamically!"
							  << std::endl;
						std::cerr << "  Pair Bond calc. for outer event from inside: " << i1 << "-" << i2 << " alpha = " 
							  << alpha << "  term = " << term
							  << "  discr = " << discr << "  part1 (neg) = " << part1 << "  part2 (pos) = " << part2
							  << "  coll. time = " << collisionTime << " time = " << timeNow << std::endl;
						std::cerr << "  Printing previous event:" << std::endl;
						evSchedule.print();
							
						errorFlag = true;
#ifdef USE_OLDQ
						ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
						hybridQ->ScheduleEvent(1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
						
					}
#ifdef USE_OLDQ
					ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j);
#else
					hybridQ->ScheduleEvent(i1, i2, timeNow+collisionTime, 1.0,  BOND_COLLISION, j, interactionNumber[i1], interactionNumber[i2]);
#endif
				}
			}
		} // end discr>0
	} // end of index matching
  } // end sum over j
}


void ProteinChain::calculateAllCollisions(int na, int nb){

  if (nb < 0){
    for (int i = na+1;i < numAtoms;i++){
      calculateCollisionPair(na, i);
    }
  } else {
    for (int i = 0;i < numAtoms;i++){
      if ( i == na ) continue;
      calculateCollisionPair(na,i); // includes na-nb pair
      if ( i == nb ) continue;
      calculateCollisionPair(nb, i);
  
    }
  } 

}

void ProteinChain::calculateCollisionPair(int i1, int i2){

  bool boundFlag = false; // set to true if in an infinite SW
  double minBondLengthsq = minBondLengthSq_hc; // hard core repulsion between non-connected pairs
  double maxBondLengthsq = INFINITY;

  //Vector r12 = pos[i2] - pos[i1] + (timeNow - localTime[i2])*vel[i2] - (timeNow - localTime[i1])*vel[i1];

  double dt1 = timeNow - localTime[i1];
  double dt2 = timeNow - localTime[i2];
  Vector r1 = pos[i1] + vel[i1]*dt1;
  Vector r2 = pos[i2] + vel[i2]*dt2;
  Vector r12 = r2 - r1;

  double r12sq = r12.nrm2();

  Vector v12 = vel[i2] - vel[i1];  
  double v12sq = v12.nrm2();
  
  double alpha= (r12*v12) /v12sq; // dot product

  double minCollisionTime = INFINITY;
  double maxCollisionTime = INFINITY;

  if ( abs(i2-i1) == 1 ) {
    // nearest neighbors
	  minBondLengthsq = minBondLengthSq1;
	  maxBondLengthsq = maxBondLengthSq1;
	  boundFlag = true;
  } else if ( abs(i2-i1) == 2 ){
	  minBondLengthsq = minBondLengthSq2;
	  maxBondLengthsq = maxBondLengthSq2;
	  boundFlag = true;
  } else {
	  minBondLengthsq =  minBondLengthSq_hc;
	  if (alpha > 0.0) {
             if (r12sq < minBondLengthsq){
                if ( ( (i1==evSchedule.idA) and (i2==evSchedule.idB) ) or ( (i1==evSchedule.idB) and (i2==evSchedule.idB) ) ) return;


                std::cerr << "  Hard core problem for indices " << i1 << "-" << i2 
                           << "  r12sq = " << r12sq << " is less than " << minBondLengthsq
                           << " diff = " << r12sq - minBondLengthsq << std::endl;
                evSchedule.print();
                //exit(1);
             }
             return; // no collision at hard sphere distance since particles moving away from one another
	 }
  }


  if ( alpha < 0.0){
    // Moving in towards minimum distance First look for a repulsion at short distance

    double term = (r12sq-minBondLengthsq)/v12sq;
    if (term < 0){
        std::cerr << " May have missed an inner bond collision since r12sq < " << minBondLengthsq << "  diff = " << r12sq - minBondLengthsq << std::endl;
        std::cerr << "  Colliding pair are " << i1 << "-" << i2 << "  v12sq = " << v12sq << std::endl;
        //exit(1);
    }

    double discr = alpha*alpha - term;
    if (discr < 0.0){
      //
      //   No repulsion since no real roots so must check for a stretching collision
      //
      if (!boundFlag) return; // no stretching collision if pair not bound in SW

      term = (r12sq-maxBondLengthsq)/v12sq; // this must be negative or zero
      discr = alpha*alpha - term; // new discriminant - this must be positive since must be inside well
      if (discr < 0.0){
	  std::cerr << "Error!  Is pair " << i1 << "-" << i2
		<< " out of bonding well?   r12 = " << sqrt(r12sq) << "  term = "
		<< term <<  "  time = " << timeNow << std::endl;
          errorFlag = true;
#ifdef USE_OLDQ
	  ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1);
#else
	  hybridQ->ScheduleEvent(-1, -1, timeNow, 0.0,  STOP, -1, -1, -1);
#endif
	  return;
      }

      double part1 = -alpha - sqrt(discr); // expect this root to be negative or zero
      double part2 = -alpha + sqrt(discr); // should be collision time to hit outer boundary

      collisionType = MAX_COLLISION;
      maxCollisionTime=part2;

#ifdef PRINT_SW
      if (timeNow > debugTime) { 
	  if ( (i1 == debugId1) or (i1 == debugId2)) {
		  std::cout << "   Indices " << i1 << "-" << i2 << " part1 neg. or zero, part2 positive." << std::endl;
		  std::cout << "  In max collision, Alpha = " << alpha << "  part1 = " << part1 << "  part2 = " 
			<< part2 << std::endl << "      term = " << term << "  maxCollision time = "
			<< maxCollisionTime << "  distance to maxBond = " 
			<< r12sq-maxBondLengthsq << " r12 = " << sqrt(r12sq) << endl; 
	  }
      }    
#endif

    } else {
      // 2 real roots exist since original discr >= 0.0
      double part1 = -alpha - sqrt(discr);
      double part2 = -alpha + sqrt(discr);
      
      //   After collision at short distance, would have alpha > 0.0, so both roots should be positive
      collisionType = MIN_COLLISION;
      minCollisionTime = part1;


#ifdef PRINT_SW
      if (timeNow > debugTime) { 
	  if ( (i1 == debugId1) or (i1 == debugId2)) {
		  std::cout << "   Indices " << i1 << "-" << i2 << " both part1 and part2 positive" << std::endl;
		  std::cout << "  In short collision, Alpha = " << alpha << "  part1 = " << part1 << "  part2 = " 
			<< part2 << std::endl << "      term = " << term << "  minCollision time = "
			<< minCollisionTime << "  distance to minBond = " 
			<< r12sq-minBondLengthsq << " r12sq = " << r12sq << endl; 
	  }
      }
      
#endif
    }
  } else {
    //
    //   Alpha is positive, and always increases: relative distance initially increasing
    //   Only look for a stretching collision    
    double term = (r12sq-maxBondLengthsq)/v12sq; // this must be negative and non-zero (alpha > 0) means not at term=0
    
    double discr = alpha*alpha-term;
    if (discr < 0.0) return; // no minimum collision possible since alpha increasing
 
    double part1 = -alpha - sqrt(discr); // should be negative and non-zero
    double part2 = -alpha + sqrt(discr);

    collisionType = MAX_COLLISION;
    maxCollisionTime = part2;


#ifdef PRINT_SW
    if (timeNow > debugTime) { 
	if ( (i1 == debugId1) or (i1 == debugId2)) {
 		std::cout << "   Indices " << i1 << "-" << i2 << "  part1 neg. and not zero, part2 positive" << std::endl;
		std::cout << "  In max collision, Alpha = " << alpha << "  part1 = " << part1 << "  part2 = "
			  << part2<< std::endl  << "      term = " << term
			  << "  maxCollision time = " << maxCollisionTime << "  distance to maxBond = " 
			  << r12sq-maxBondLengthsq << " r12 = " << sqrt(r12sq) << endl; 
	}
    }    
#endif
    
  }

  
  if (maxCollisionTime < INFINITY) {
#ifdef USE_OLDQ
	  ScheduleEvent(i1, i2, timeNow+maxCollisionTime, 0.0, MAX_COLLISION, -1);
#else
	  hybridQ->ScheduleEvent(i1, i2, timeNow+maxCollisionTime, 0.0, MAX_COLLISION, -1, interactionNumber[i1], interactionNumber[i2]);
#endif

#ifdef PRINT
	  if ( timeNow > debugTime ){
		std::cout << std::setprecision(12) 
			<< " Scheduling:  " << i1 << "-" << i2
			<< "  max coll. time = " 
			<< timeNow + maxCollisionTime << " r12sq = " << r12sq << std::endl;
                checkCollision(i1, i2, timeNow+maxCollisionTime, maxBondLengthsq);
          }
#endif

  } else if (minCollisionTime < INFINITY) {
#ifdef USE_OLDQ
	  ScheduleEvent(i1, i2, timeNow+minCollisionTime, 0.0,  MIN_COLLISION, -1); 
#else
	  hybridQ->ScheduleEvent(i1, i2, timeNow+minCollisionTime, 0.0,  MIN_COLLISION, -1, interactionNumber[i1], interactionNumber[i2]);
#endif

#ifdef PRINT
	  if ( timeNow > debugTime ) {   
		  std::cout << std::setprecision(12) << " Scheduling:  " << i1 << "-" << i2
				<< "  min coll. time = " << timeNow + minCollisionTime << " r12sq = " << r12sq << std::endl;
		 checkCollision(i1, i2, timeNow+minCollisionTime, minBondLengthsq);
          }
#endif
    
  }

}

void ProteinChain::checkCollision(int collider1, int collider2, double cTime, double dsq){
    double dt1 = cTime - localTime[collider1];
    double dt2 = cTime - localTime[collider2];
   
    Vector r1 = pos[collider1] + vel[collider1]*dt1; 
    Vector r2 = pos[collider2] + vel[collider2]*dt2;

    Vector r12 = r2 - r1;
    double rsq = r12.nrm2();
    double diff = rsq - dsq;
    std::cout << "  In checkCollision between " << collider1 << "-" << collider2 <<  "  diff = " 
              << diff << " cTime = " << cTime << "  dt1 = " << dt1 << "   dt2 = " << dt2 << "  r12sq = " << rsq 
              << std::endl << "           r1 = " << r1 << "  r2 = " << r2 << std::endl;
  

}
