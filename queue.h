#ifndef QUEUE_H
#define QUEUE_H
#include "event.h"
#include <iostream>
#include <queue>
#include <vector>
#include <random>
#include <memory>

const int nlists = 101;
const int poolSize = 50000;

struct compareEvents
{
    bool operator()(const std::shared_ptr<Event> lhs, const std::shared_ptr<Event> rhs)
    {
        return lhs->time > rhs->time;
    }
};



template <class T, class S, class C>
S& Container(std::priority_queue<T, S, C>& q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static S& Container(std::priority_queue<T, S, C>& q) {
            return q.*&HackedQueue::c;
        }
    };
    return HackedQueue::Container(q);
}

class HybridQueue 
{
    private:
        std::default_random_engine generator;
        std::vector<std::shared_ptr<Event>> linearList[nlists+1]; // will hold events not in priority queue, time slice of events

        int currentIndex; // reference to linear list index to process
        int baseIndex; // used to go through linear list in periodic fashion
        int numInList; // number of events put into linear list
        double listIntervalTime;
        int listUpdateCount;
        double averageListSize;
        double scale; // determines balance between priority queue and linear list

        std::shared_ptr<Event> currentEvent; // allow non-friend access to current event?

        std::vector<std::shared_ptr<Event>> eventPool; // container for pool of event pointers to use

        std::priority_queue<std::shared_ptr<Event>, std::vector<std::shared_ptr<Event>>, compareEvents> hQ; // name of actual ordered queue

        std::shared_ptr<Event> drawFromPool(); // Not needed if no pool (dynamic allocation of events)

    public:
        HybridQueue(double scale = 500);
        ~HybridQueue();

        void printQueue();
        void printTopQueue();
        
        void resetEvent(std::shared_ptr<Event> ev);
        void reset();
        void processOverflow(); // checks to see if distant events are still valid, and feeds into linear list
        void drawQueue(); // draws next event ptr and places in currentEvent
        void drawQueue(Event& ev, std::vector<int>& interactionNumber);
        void replenishQueue();
        void replenishQueue(std::vector<int>& interactionNumber);
        void insertInEventQueue(std::shared_ptr<Event> ev);
        void insertInEventQueue(std::shared_ptr<Event> ev, std::vector<int> &interactionNumber);

        void ScheduleEvent(int idA, int idB, double tEvent,
                            double deltaU, int collType, int bondIndex, int interactionA, int interactionB);
        void ScheduleEvent(const Event& event);

        void runTest(int numInitialEvents, int numEvents,double maxTime);

        double getAverageQSize(){ return double(averageListSize)/double(listUpdateCount); }
        


};
#endif
