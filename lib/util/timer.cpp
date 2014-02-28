#include <excafe/util/timer.hpp>
#include <time.h>
#include <cstddef>

namespace excafe
{

namespace util
{

Timer::Timer() : running(false), seconds(0), nanoseconds(0)
{
}

void Timer::start()
{
  if (!running)
  {
    running = true;
    clock_gettime(CLOCK_MONOTONIC, &begin);
  }
}

void Timer::stop()
{
  if (running)
  {
    running = false;

    timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);

    seconds += end.tv_sec - begin.tv_sec;
    nanoseconds += end.tv_nsec - begin.tv_nsec;

    adjustNanoseconds();
  }
}

void Timer::adjustNanoseconds()
{
  while(nanoseconds < 0)
  {
    nanoseconds += billion;
    --seconds;
  }

  while(nanoseconds >= billion)
  {
    nanoseconds -= billion;
    ++seconds;
  }
}

double Timer::getSeconds() const
{
  return static_cast<double>(seconds) + static_cast<double>(nanoseconds)/billion;
}

long Timer::getIntegerSeconds() const
{
  return seconds;
}

double Timer::getFractionalSeconds() const
{
  return static_cast<double>(nanoseconds)/billion;
}

}

}
