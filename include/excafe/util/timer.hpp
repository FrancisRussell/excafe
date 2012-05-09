#ifndef EXCAFE_UTIL_TIMER_HPP
#define EXCAFE_UTIL_TIMER_HPP

#include <sys/time.h>

namespace excafe
{

namespace util
{

class Timer
{
private:
  static const long million = 1000000;

  bool running;
  timeval begin;
  long seconds;
  long microseconds;

  void adjustMicroseconds();

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
