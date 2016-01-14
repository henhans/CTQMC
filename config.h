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

    // operators
    // copy configurations
    Time_config& operator=( Time_config& rhs);

    // print out time configuration for each flavor
    void print_config(int fl_);
    // check if the time configuration for flavor i is empty
    bool is_empty(int fl_) { return ( t_start[fl_].empty() && t_end[fl_].empty() ); };
    // get number of flavor
    int get_flavor() { return flavor; };
    // get current order for flavor
    int get_pertur_order(int fl_) { return t_start[fl_].size(); };
    // get t_start element
    double get_t_start_at(int fl_, int indx) { return t_start[fl_][indx]; };
    // get t_end element
    double get_t_end_at(int fl_, int indx) { return t_end[fl_][indx]; } ;

    // insert time sector
    pair<bool, int> try_insert_start_time ( int fl_, double ts_);//check if the propose starting time can be inserted
    //then return the bolean constant if ts ca be insert and the index for insertion
    void insert_start_time ( int fl_, int index_ , double ts_);//insert the starting time
    pair<double, int> find_next_start_time ( int fl_, int index_);//find the next starting time for end time insertion
    //return the next starting time for generate random te and the index for insertion
    void insert_end_time ( int fl_, int index_ ,double te_);//insert the end time

    // remove time sector
    void remove_time_sector(int fl_, int index_);

    
};

Time_config& Time_config::operator = ( Time_config& rhs ) {

    flavor = rhs.get_flavor();

    for(int ifl=0; ifl<flavor ; ifl++) {// for each flavor
        t_start[ifl].clear();// clear t_start
        t_end[ifl].clear();// clear t_end
        for(int i=0; i<rhs.get_pertur_order(ifl); i++) {
            t_start[ifl].push_back(rhs.t_start[ifl][i]);// copy t_start
            t_end[ifl].push_back(rhs.t_end[ifl][i]);// copy t_end
        }
    }

    return *this;
};


void Time_config::print_config( int fl_) {
    clog << "t start config:" << endl;
    for (deque<double>::iterator it = t_start[fl_].begin(); it!=t_start[fl_].end(); ++it)
        clog << ' ' << setprecision(9) << *it;
    cout << '\n';

    clog << "t end config:" << endl;
    for (deque<double>::iterator it = t_end[fl_].begin(); it!=t_end[fl_].end(); ++it)
        clog << ' ' << setprecision(9) << *it;
    cout << '\n';

}

pair<bool, int> Time_config::try_insert_start_time( int fl_, double ts_) {
    if( t_start[fl_].empty() ) return make_pair(true, -1);// return -1 if the deque is empty
    else{
        //bool accept = false;
        //int index = -1;

        //find if the ts insertion is allowed and return the acceptance condition and the index
        int ts_size = t_start[fl_].size();
        //check if the propose time sector is locate at the periodic sector
        if( (t_end[fl_][ts_size-1] < t_start[fl_][ts_size-1]) && (ts_ > t_start[fl_][ts_size-1] || ts_ <
            t_end[fl_][ts_size-1]) )
            return make_pair(false,-1);
        //clog << "here?" <<endl;

        for(int i=0; i < ts_size; i++) {// then locate the insertion time
            //clog <<i<<" "<< ts_size<< " " <<ts_ <<" "<< t_start[fl_][i] << " "<< t_end[fl_][(i-1+ts_size)%ts_size] << endl;
            if( ts_ == t_start[fl_][i] )
                return make_pair(false,-1);
            else if (ts_< t_start[fl_][i] && i>0 && ts_< t_end[fl_][(i-1)]/* && i>0*/) {
                //clog << "here1?" << endl;
                return make_pair(false,-1);
            }
            else if(ts_ < t_start[fl_][i] && ts_ > t_end[fl_][(i-1+ts_size)%ts_size] ) {
                //clog << "here2?" << endl;
                return make_pair(true,i);
            }
            else if (ts_ < t_start[fl_][i] && ts_ < t_end[fl_][(i-1+ts_size)%ts_size] ) {
                //clog << "here3?" << endl;
                return make_pair(true,i);
            }
            //clog << "here4?" << endl;
        }
        //clog << "here?" << endl;
        if(ts_ > t_end[fl_][ts_size-1]) return make_pair(true,ts_size);
        else return make_pair(false,-1);
    };
};

void Time_config::insert_start_time ( int fl_, int index_ ,double ts_ ) {
    if( index_==-1 ) t_start[fl_].push_front(ts_);
    else if( index_== t_start[fl_].size() ) t_start[fl_].push_back(ts_);
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

void Time_config::remove_time_sector(int fl_, int index_) {
    t_start[fl_].erase( t_start[fl_].begin() + index_ );
    t_end[fl_].erase( t_end[fl_].begin() + index_ );
};

