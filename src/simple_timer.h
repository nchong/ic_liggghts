#include <sys/time.h>
#include <unistd.h>

class SimpleTimer {
  public:
    SimpleTimer() : _total_time(0.0f) {}
    ~SimpleTimer() {}

    inline void start() { gettimeofday(&_start, NULL); };
    inline void stop() { gettimeofday(&_end, NULL); };
    inline void add_to_total() { 
      long seconds = _end.tv_sec - _start.tv_sec;
      long useconds = _end.tv_usec - _start.tv_usec;
      _total_time += ((seconds) * 1000 + useconds/1000.0) + 0.5;
    };
    inline double total_time() { return _total_time; };
    inline void reset() { _total_time = 0.0; };

  private:
    struct timeval _start, _end;
    double _total_time;
};
