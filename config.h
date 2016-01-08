#include <deque>
using namespace std;

typedef map<int, deque<double> > map_of_deque; // data structure for a map of queue

// class for time configuration and its operators(add time, remove time...)
class Time_config{

    int flavor; // number of flavor
    map_of_deque t_start;// time start configuration: key is flavor, element is queue of configuration for each flavor.
    map_of_deque t_end;// time end configuration:...

    public:
    // constructor
    Time_config(int fl_){
        flavor = fl_;
        for (int i=0; i<fl_; i++) {
            clog << "initial empty time queue for flavor:" << i << endl;
            deque<double> empty;// temperory empty queue
            t_start[i] = empty;
            t_end[i] = empty;
            //clog << "size of queue is " << t_config[i].size() << endl;
        };
    };

    // print out time configuration for each flavor
    void print_config(int fl_);
    // check if the time configuration for flavor i is empty
    bool is_empty(int fl_) { return ( t_start[fl_].empty() && t_end[fl_].empty() ); };
    // get number of flavor
    int get_flavor() { return flavor; };

    // insert time
    pair<bool, int> try_insert_start_time ( int fl_, double ts_);//check if the propose starting time can be inserted
    //then return the bolean constant if ts ca be insert and the index for insertion
    void insert_start_time ( int fl_, int index_ , double ts_);//insert the starting time
    pair<double, int> find_next_start_time ( int fl_, int index_);//find the next starting time for end time insertion
    //return the next starting time for generate random te and the index for insertion
    void insert_end_time ( int fl_, int index_ ,double te_);//insert the end time

    // remove time

};

void Time_config::print_config( int fl_) {
    clog << "t start config:" << endl;
    for (deque<double>::iterator it = t_start[fl_].begin(); it!=t_start[fl_].end(); ++it)
        clog << ' ' << *it;
    cout << '\n';

    clog << "t end config:" << endl;
    for (deque<double>::iterator it = t_end[fl_].begin(); it!=t_end[fl_].end(); ++it)
        clog << ' ' << *it;
    cout << '\n';

}

pair<bool, int> Time_config::try_insert_start_time( int fl_, double ts_) {
    if( t_start[fl_].empty() ) return make_pair(true, -1);// return -1 if the deque is empty
    else{
        //bool accept = false;
        //int index = -1;
        //find if the ts insertion is allowed and return the acceptance condition and the index
        int ts_size = t_start[fl_].size();
        for(int i=0; i < ts_size; i++) {
            //clog << t_start[fl_][i] << " "<< t_end[fl_][(i-1+ts_size)%ts_size] << endl;
            if(ts_ < t_start[fl_][i] && ts_ > t_end[fl_][(i-1+ts_size)%ts_size] )
                return make_pair(true,i);
            else if (ts_ < t_start[fl_][i] && ts_ < t_end[fl_][(i-1+ts_size)%ts_size] )
                return make_pair(true,i);
        }
        if(ts_ > t_end[fl_][ts_size-1]) return make_pair(true,ts_size-1);
        else return make_pair(false,-1);
    };
};

void Time_config::insert_start_time ( int fl_, int index_ ,double ts_ ) {
    if( index_==-1 ) t_start[fl_].push_front(ts_);
    else if( index_== t_start[fl_].size()-1 ) t_start[fl_].push_back(ts_);
    else{// insert ts at specific position
        deque <double>::iterator Iter=t_start[fl_].begin();
        for(int i=0; i<index_; i++) Iter++;
        t_start[fl_].insert(Iter, ts_);
    };
};

pair<double, int> Time_config::find_next_start_time ( int fl_, int index_) {
    if( index_ ==-1 ) return make_pair(t_start[fl_][index_+1], -1);// return next time and index -1
    else{
        int ts_size = t_start[fl_].size();
        return make_pair( t_start[fl_][(index_+1)%ts_size], (index_+1)%ts_size);
    };
};

void Time_config::insert_end_time( int fl_, int index_, double te_) {
    if(index_ == -1 ) t_end[fl_].push_front(te_);
    else{
        deque <double>::iterator Iter=t_end[fl_].begin();
        for(int i=0; i<index_; i++) Iter++;
        t_end[fl_].insert(Iter, te_);
    };

}
