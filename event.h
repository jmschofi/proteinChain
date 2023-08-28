#ifndef EVENT_H
#define EVENT_H
#include <limits>
#include <numeric>
#include <iostream>
#include <iomanip>

template <typename T>
    struct counter
    {
        counter()
        {
            objects_created++;
            objects_alive++;
        }

        counter(const counter&)
        {
             objects_created++;
             objects_alive++;
        }   

    //protected:
        virtual ~counter()
        {
            --objects_alive;
        }
        static int objects_created;
        static int objects_alive;
    };

template <typename T> int counter<T>::objects_created( 0 );
template <typename T> int counter<T>::objects_alive( 0 );

class Event : counter<Event> {

  public:
    int index;
    double time; // will hold key to sorting events   
    double deltaU;
    int idA, idB,  collType;
    int interactionA, interactionB, bondIndex;

    int left, right, up, circAR; // old queue compatibility
    int qIndex, next, previous; // used for hybrid tree

    void print(){
      int collision1 = idA;
      int collision2 = idB;
      if (collision1 > collision2) {
                  collision1 = idB;
                  collision2 = idA;
      }
      std::cout <<  " Pair " << collision1 << "-" << collision2 << "  collType is " << collType
                                << "  with dU = " << deltaU << " bond index = " << bondIndex <<  " at time "
                                << std::setprecision(8) << std::setw(8) << time << "  interactonA = " << interactionA
                                << " interactionB = " << interactionB << std::endl;
    }

    Event() : idA(-1), idB(-1), interactionA(-1), interactionB(-1), deltaU(0.0), time(std::numeric_limits<double>::max() ) {}

    Event(int id, double t=std::numeric_limits<double>::max()) : deltaU(0.0), idA(-1), idB(-1), interactionA(-1), interactionB(-1),  
                  index(id), time(t) {}
    ~Event() { }

    void printEvent(){ std::cout << " Event index " << index << " is at time " 
            << std::setprecision(16) << time << "  idA = " << idA << " idB = " << idB 
            << "  collType = " << collType << " interactionA = " << interactionA << " interactionB = " << interactionB << std::endl; }

};

#endif
