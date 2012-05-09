#include <excafe/util/timer.hpp>
#include <sys/time.h>
#include <cstddef>

namespace excafe
{

namespace util
{

Timer::Timer() : running(false), seconds(0), microseconds(0)
{
}

void Timer::start()
{
  if (!running)
  {
    running = true;
    gettimeofday(&begin, NULL);
  }
}

void Timer::stop()
{
  if (running)
  {
    running = false;

    timeval end;
    gettimeofday(&end, NULL);

    seconds += end.tv_sec - begin.tv_sec;
    microseconds += end.tv_usec - begin.tv_usec;

    adjustMicroseconds();
  }
}

void Timer::adjustMicroseconds()
{
  while(microseconds < 0)
  {
    microseconds += million;
    --seconds;
  }

  while(microseconds >= million)
  {
    microseconds -= million;
    ++seconds;
  }
}

double Timer::getSeconds() const
{
  return static_cast<double>(seconds) + static_cast<double>(microseconds)/million;
}

long Timer::getIntegerSeconds() const
{
  return seconds;
}

double Timer::getFractionalSeconds() const
{
  return static_cast<double>(microseconds)/million;
}

}

}
