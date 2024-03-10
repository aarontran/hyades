#ifndef TIMER_H
#define TIMER_H

typedef struct timekeep {
  double t;
  double tbuf;
  double ttot;
  //double ttotal;  // not implemented
  //int ncalls;
  char name[100];  // see also TimerArray member "max_name_len"
} timekeep_t;

// Name "timer_t" is already used in C standard library

class TimerArray {

  private:
    const static int max_name_len = 100;

  public:
    TimerArray(int ntmax_);

    int ntmax;   // max number of timers
    int nt;      // current number of timers
    timekeep_t* t0;

    // Data

    // Low-level subroutines
    double wallclock();
    timekeep_t* seek(const char* name);

    // High-level subroutines
    void add(const char* name);
    void tic(const char* name);  // not implemented
    void toc(const char* name);
    double read(const char* name);
    double read_total(const char* name);
    double flush(const char* name);
    void flush_all();

};

#endif  // TIMER_H
