#include <assert.h>
#include <stdlib.h>  // for malloc, free
#include <string.h>
#include <sys/time.h>

#include "timer.h"

// Constructor
TimerArray::TimerArray(int ntmax_
):
  ntmax(ntmax_)
{
  nt = 0;
  t0 = (timekeep_t*) malloc( ntmax_*sizeof(timekeep_t) );
}

timekeep_t* TimerArray::seek(const char* name) {
  timekeep_t* watch = t0;
  int ii = 0;
  while (ii < nt) {
    if (strcmp(name, watch->name) == 0) break;
    ++ii;
    ++watch;
  }
  assert(ii < nt);  // error out if cannot find timer
  return watch;
}

void TimerArray::add(const char* name) {
  assert(nt < ntmax);
  timekeep_t* watch = &(t0[nt]);
  watch->t    = -1;  // IMPORTANT sentinel value, < 0 means not ticking.
  watch->tbuf = 0;
  watch->ttot = 0;
  strncpy(watch->name, name, max_name_len);
  nt += 1;
}

// cannot already be ticking (don't tic + tic without intervening toc)
void TimerArray::tic(const char* name) {
  timekeep_t* watch = seek(name);
  assert(watch->t < 0);
  watch->t = wallclock();
}

// "toc" must be preceded by "tic"
void TimerArray::toc(const char* name) {
  timekeep_t* watch = seek(name);
  assert(watch->t > 0);
  watch->tbuf += (wallclock() - watch->t);
  watch->t = -1;  // reset for next "tic"
}

// "read" must be preceded by "tic"
double TimerArray::read(const char* name) {
  timekeep_t* watch = seek(name);
  assert(watch->t > 0);
  return wallclock() - watch->t;
}

// must come after "toc" and "flush"
double TimerArray::read_total(const char* name) {
  timekeep_t* watch = seek(name);
  assert(watch->t < 0);  // cannot be ticking
  assert(watch->tbuf == 0);  // must be flushed
  return watch->ttot;
}

// "flush" must be preceded by "toc" to guard against mis-use
// OK to call this twice in a row (second flush is a no-op)
double TimerArray::flush(const char* name) {
  timekeep_t* watch = seek(name);
  assert(watch->t < 0);
  double result = watch->tbuf;
  watch->tbuf = 0;
  watch->ttot += result;
  return result;
}

// "flush" must be preceded by "toc" to guard against mis-use
// OK to call this twice in a row (second flush is a no-op)
void TimerArray::flush_all() {
  timekeep_t* watch = t0;
  int ii = 0;
  while (ii < nt) {
    assert(watch->t < 0);
    double result = watch->tbuf;
    watch->tbuf = 0;
    watch->ttot += result;
    ++ii;
    ++watch;
  }
  return;
}

// TODO use more modern timer functions --ATr,2024feb25
// https://people.cs.rutgers.edu/~pxk/416/notes/c-tutorials/gettime.html
double TimerArray::wallclock(void) {
  //https://stackoverflow.com/questions/5362577/c-gettimeofday-for-computing-time
  //https://people.cs.rutgers.edu/~pxk/416/notes/c-tutorials/gettime.html
  struct timeval tv[1];
  gettimeofday( tv, NULL );
  return (double)(tv->tv_sec) + 1e-6*(double)(tv->tv_usec);
}
