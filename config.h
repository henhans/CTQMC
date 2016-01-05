#include <queue>
using namespace std;

typedef map<int, queue<double> > map_of_queue; // data structure for a map of queue

// class for time configuration and its operators(add time, remove time...)
class Time_config{

   int flavor; // number of flavor
   map_of_queue t_start;// time start configuration: key is flavor, element is queue of configuration for each flavor.
   map_of_queue t_end;// time end configuration:...

public:
   // constructor
   Time_config(int fl_){
     flavor = fl_;
     for (int i=0; i<fl_; i++) {
       clog << "initial empty time queue for flavor:" << i << endl;
       queue<double> empty;// temperory empty queue
       t_start[i] = empty;
       t_end[i] = empty;
       //clog << "size of queue is " << t_config[i].size() << endl;
     };
   };

   // get number of flavor
   int get_flavor() {
     return flavor;
   }
   // insert time

   // remove time
};
