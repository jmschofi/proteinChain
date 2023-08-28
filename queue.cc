#include "queue.h"

//   See Gerald Paul, J. Comp. Phys. 221, 615 (2007)

//#define VERBOSE
//#define PRINT_AVERAGE_QSIZE

HybridQueue::HybridQueue(double scale) : currentIndex(0), baseIndex(0), eventPool(0), currentEvent(0),
                        listIntervalTime(1), averageListSize(0), listUpdateCount(0), scale(scale) {

    for (int i=0;i<=nlists;i++) linearList[i].reserve(50);
    //
    //  create a pool of event pointers for use to avoid copy constructors
    //
    for (int i=0;i<poolSize;i++){
        std::shared_ptr<Event> ev_ptr(new Event(i));
        eventPool.push_back(ev_ptr);
    }

}

HybridQueue::~HybridQueue(){
    std::cout << " Destructor for Hybrid Queue called." << std::endl;
    //
    //  Clean up the event Pool
    std::cout << "alive:  Events: " << counter<Event>::objects_alive << std::endl;
}

std::shared_ptr<Event> HybridQueue::drawFromPool(){
    if (eventPool.empty()){
        std::cout << " !!!!!!!!  %%%%%%   ######  Expanding event pool:" << std::endl;
        for (int i=0;i<poolSize;i++){
            std::shared_ptr<Event> ev_ptr(new Event(i) );
            eventPool.push_back(ev_ptr);
        }
    }
    std::shared_ptr<Event> ev = eventPool.back();
    eventPool.pop_back();
    return ev;
}

void HybridQueue::runTest(int numInitialEvents, int numEvents, double maxTime){ 
    //
    listIntervalTime = maxTime/double(nlists-1);

    //  Randomly create event times
    std::uniform_real_distribution<double> distribution(0.0, maxTime/10.0);
    // create events to analyze and place in queue or in linear list
    int count = 0;
    while (count < numInitialEvents){
        double ltime = distribution(generator);
        std::shared_ptr<Event> ev = drawFromPool();
        ev->index = count;
        ev->time = ltime;

        insertInEventQueue(ev);
        count++;
    }

    //
    //  Print out initial queue
    //
#ifdef VERBOSE
    printQueue();
#endif
   
    //
    //  Now draw events
    double prevTime = 0.0;
    for (int i=0;i<numEvents;i++){
        drawQueue();

#ifdef VERBOSE
        std::cout << "  Event " << i << " drew event:" << std::endl;
        currentEvent->printEvent();
#endif

        if (currentEvent->time < prevTime){
            std::cout << "  Error here!  previous time = " << prevTime << "  currentTime = " 
                            << currentEvent->time << std::endl;
            exit(1);
        }
        prevTime = currentEvent->time;
        // 
        // Event processed here: nothing done for now
        //
        eventPool.push_back( currentEvent );  // return point to pool

        //  Replace with new event
        count++;
        double ltime = distribution(generator);
        std::shared_ptr<Event> newEvent = drawFromPool();
        newEvent->index = count;
        newEvent->time = prevTime + ltime; 

        insertInEventQueue(newEvent);
    }

    std::cout << "  average size of list is " << averageListSize/double(listUpdateCount) << std::endl;
    
    std::cout << "created: "
             << " Events: " << counter<Event>::objects_created
             << std::endl;
    std::cout << "alive: "
             << " Events: " << counter<Event>::objects_alive
             << std::endl;
}

void HybridQueue::drawQueue(){
    if (hQ.empty() ) replenishQueue();

    currentEvent = hQ.top(); 
    hQ.pop(); // remove event from queue

}


void HybridQueue::drawQueue(Event& ev, std::vector<int>& interactionNumber){

    bool valid = false;
    while ( !valid ){    
        if (hQ.empty() ) replenishQueue(interactionNumber); // allow filtering

        currentEvent = hQ.top(); 
        hQ.pop(); // remove event from queue

        int idA = currentEvent->idA;
        int idB = currentEvent->idB;
        if (idA > -1){
            if (currentEvent->interactionA != interactionNumber[idA] || currentEvent->interactionB != interactionNumber[idB]) {
                valid = false;
                resetEvent( currentEvent ); // reset all values to default
                eventPool.push_back( currentEvent );  // return point to pool
            } else {
                valid = true;
            }
        } else {
            valid = true;
        }

    }

    ev.idA = currentEvent->idA;
    ev.idB = currentEvent->idB;
    ev.time = currentEvent->time;
    ev.deltaU = currentEvent->deltaU;
    ev.collType = currentEvent->collType;
    ev.bondIndex = currentEvent->bondIndex;
    ev.interactionA = currentEvent->interactionA;
    ev.interactionB = currentEvent->interactionB;

    eventPool.push_back( currentEvent ); // return stored event to pool
#ifdef VERBOSE
    std::cout << " In new drawQueue, event drawn is:" << std::endl;
    ev.print();
#endif

}

void HybridQueue::replenishQueue(){
    while (hQ.empty() ){
        //
        //   Tree is empty
        //
        currentIndex++; // move to next linked list
        if (currentIndex == nlists){
            currentIndex = 0;
            baseIndex += nlists;
            processOverflow();
        }
        //
        //   Feed list into tree
        //
#ifdef VERBOSE
        std::cout << "  *****   Replenishing q with " << linearList[currentIndex].size() 
              << " elements in list " << currentIndex << std::endl;
#endif

        listUpdateCount++;
        averageListSize += linearList[currentIndex].size();
        
        while ( !linearList[currentIndex].empty() ) {
            hQ.push( linearList[currentIndex].back() );
            linearList[currentIndex].pop_back();
        }   
        
    }
}

void HybridQueue::replenishQueue(std::vector<int> & interactionNumber){

#ifdef VERBOSE
    std::cout << std::endl << " Replenish Q called. Current state of lists is:" << std::endl;
    printQueue();
#endif

    while (hQ.empty() ){
        //
        //   Tree is empty
        //
        currentIndex++; // move to next linked list
        if (currentIndex == nlists){
            currentIndex = 0;
            baseIndex += nlists;
            processOverflow();
        }
        //
        //   Feed list into tree
        //
#ifdef VERBOSE
        std::cout << "  *****   Replenishing q with " << linearList[currentIndex].size() 
              << " elements in list " << currentIndex << std::endl;
#endif
        
        while ( !linearList[currentIndex].empty() ) {
            //
            //  Filter for validity
            //
            std::shared_ptr<Event> ev = linearList[currentIndex].back();
            if ( ev->idA > -1){
                int idA = ev->idA;
                int idB = ev->idB;
                if (ev->interactionA != interactionNumber[idA] || ev->interactionB != interactionNumber[idB]) {
                    // no longer valid so return to ppol
#ifdef VERBOSE
                    std::cout << "  Event is no longer valid:" << std::endl;
                    ev->print();
#endif
                    linearList[currentIndex].pop_back();
                    eventPool.push_back( ev ); // return stored event to pool
                } else {
                    hQ.push( linearList[currentIndex].back() );
                    linearList[currentIndex].pop_back();
                }
            } else {
                    hQ.push( linearList[currentIndex].back() );
                    linearList[currentIndex].pop_back();
            }

        }   
        
    }

    listUpdateCount++;
    averageListSize += hQ.size();

    
#ifdef VERBOSE
    std::cout << " average Q size is " << double(averageListSize)/double(listUpdateCount) << std::endl;
    std::cout << "After insertion, top queue is:" << std::endl;
    printTopQueue();
#endif

#ifdef PRINT_AVERAGE_QSIZE
    std::cout << " average Q size is " << double(averageListSize)/double(listUpdateCount) << std::endl;
#endif
}

void HybridQueue::insertInEventQueue(std::shared_ptr<Event> ev){
    int insertI = int( scale*ev->time - baseIndex);
    if (insertI > (nlists-1) ){
        //
        //  account for wrap around
        //
        insertI -= nlists; 

        if (insertI >= currentIndex-1) {
            insertI = nlists; // catch potential overflow and put at end of list
#ifdef VERBOSE
            std::cout << "  Inserting event ";
            ev->printEvent();
            std::cout << " into overflow." << std::endl;
#endif
        }

    }
#ifdef VERBOSE
    std::cout << "  New event of index " << ev->index << " at time " << ev->time 
               << " inserted at index " << insertI << " current index is " << currentIndex << std::endl;
#endif

    if (insertI == currentIndex){
        //
        //   put in Tree 
        //
        hQ.push(ev);
        //std::cout << "  Event put directly into Q!" << std::endl;
    } else { 
        //
        //  put into linear list
        //
        numInList++;
        linearList[insertI].push_back(ev);
    }
}

void HybridQueue::processOverflow(){
    //  
    //  overflow is in linear list at location nlists
    //  Check to see if event is no longer in overflow
    //
    int numChecks = linearList[nlists].size();
#ifdef VERBOSE
    std::cout << std::endl << "  *****  Process overflow called!  Current number of checks is " << numChecks << std::endl;
#endif
    std::vector<std::shared_ptr<Event>> overFlowList;
   
    // go back through overflow to decide if event needs to be taken out of overflow
    int retained = 0;
    for (int i=0; i < numChecks;i++) {
       
        std::shared_ptr<Event> ev = linearList[nlists][i]; // start at back
        int insertI = int( scale*ev->time - baseIndex) - nlists;

        if (insertI < currentIndex) {
            // insert into normal lists and queue only if not still in overflow
            // Doesn't affect current linear list since it is only for overflow
#ifdef VERBOSE
            std::cout << " Inserting event out of overflow at index " << insertI+nlists << " current index = "
                     << currentIndex << std::endl;
            ev->printEvent();
            std::cout << std::endl;
#endif
            insertInEventQueue(ev);  // will be retained if still valid in insertion process
           
        } else {
#ifdef VERBOSE
            std::cout << " Retaining event in overflow! " << std::endl;
            ev->printEvent();
#endif
            overFlowList.push_back( ev );  // retain in overflow list
            retained++;
#ifdef VERBOSE
            std::cout << "  Number retained = " << retained << std::endl;
#endif
        }
    
    }
    //  Is line below needed?
    linearList[nlists].clear();
    for (int i=0;i<overFlowList.size();i++) linearList[nlists].push_back( overFlowList[i] );

    //linearList[nlists] = overFlowList;  // should be resized
    
#ifdef VERBOSE
    std::cout << "  Retained " << retained << " = " << linearList[nlists].size() << " events in overflow. " << std::endl;
#endif

}

void HybridQueue::printTopQueue(){
    std::vector<std::shared_ptr<Event>> &tasks = Container(hQ);
    for (std::vector<std::shared_ptr<Event>>::iterator i=tasks.begin();i!=tasks.end();i++) (*i)->print();
}

void HybridQueue::printQueue(){
    std::cout << "  Current state of Queue of size " << hQ.size() << std::endl;
    std::vector<std::shared_ptr<Event>> &tasks = Container(hQ);
    for (std::vector<std::shared_ptr<Event>>::iterator i=tasks.begin();i!=tasks.end();i++) (*i)->printEvent();

    std::cout << "  Linear lists: unsorted by time:" << std::endl;
    for (int i = currentIndex; i < currentIndex + nlists; i++){

        int index = i;
        if (index > (nlists-1) ) index -= nlists;
        //if (i >= currentIndex-1) i=nlists;

        std::cout << "  Local list " << index << " contains " << linearList[index].size()
                  << " elements" << std::endl;
        for (int j=0; j < linearList[index].size(); j++) {
            linearList[index][j]->printEvent();
        }

    }   

    std::cout << " Overflow list is:" << std::endl;
    for (int j=0;j<linearList[nlists].size();j++) linearList[nlists][j]->printEvent();


}

void HybridQueue::ScheduleEvent(int idA, int idB, double tEvent,
                                 double deltaU, int collType, int bondIndex,
                                 int interactionA, int interactionB){

    std::shared_ptr<Event> ev = drawFromPool(); // get pointer from pool
    ev->idA = idA;
    ev->idB = idB;
    ev->time = tEvent;
    ev->deltaU = deltaU;
    ev->collType = collType;
    ev->bondIndex = bondIndex;
    if (idA > -1)  ev->interactionA = interactionA;
    if (idB > -1)  ev->interactionB = interactionB;

#ifdef VERBOSE
    std::cout << " In new ScheduleEvent, event scheduled was:" << std::endl;
    ev->print();
#endif

    insertInEventQueue(ev);


}

void HybridQueue::ScheduleEvent(const Event& event){

    std::shared_ptr<Event> ev = drawFromPool(); // get pointer from pool
    ev->idA = event.idA;
    ev->idB = event.idB;
    ev->time = event.time;
    ev->deltaU = event.deltaU;
    ev->collType = event.collType;
    ev->bondIndex = event.bondIndex;
    ev->interactionA = event.interactionA;
    ev->interactionB = event.interactionB;

    insertInEventQueue(ev);
}

void HybridQueue::resetEvent(std::shared_ptr<Event> ev){
    ev->idA = -1;
    ev->idB = -1;
    ev->time = std::numeric_limits<double>::max();
    ev->deltaU = 0.0;
    ev->interactionA = 0;
    ev->interactionB = 0;
    ev->bondIndex = -1;
}

void HybridQueue::reset(){

    //  Move all events back into pool
 
    for (int i=0;i<=nlists;i++){
        while ( !linearList[i].empty() ){
            std::shared_ptr<Event> ev = linearList[i].back();
            resetEvent(ev);
            eventPool.push_back(ev);
            linearList[i].pop_back();
        }
    }
    
    //  clean up priority queue
    while ( !hQ.empty() ){     
        std::shared_ptr<Event> ev = hQ.top();
        resetEvent(ev);
        eventPool.push_back( ev );  // return event to pool
        hQ.pop(); // remove event from queue
    }

    currentIndex = 0;
    baseIndex = 0;

    listUpdateCount = 0;
    averageListSize = 0.0;
    
}
