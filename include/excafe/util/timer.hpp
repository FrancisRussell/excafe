#ifndef EXCAFE_UTIL_TIMER_HPP
#define EXCAFE_UTIL_TIMER_HPP

#include <time.h>

namespace excafe
{

namespace util
{

class Timer
{
private:
  static const long billion = 1000000000l;

  bool running;
  timespec begin;
  long seconds;
  long nanoseconds;

  void adjustNanoseconds();

public:
  Timer();
  void start();
  void stop();
  double getSeconds() const;
  long getIntegerSeconds() const;
  double getFractionalSeconds() const;
};

}

}

#endif
