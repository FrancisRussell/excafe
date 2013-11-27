#ifndef EVENT_SET_HPP
#define EVENT_SET_HPP

#include <excafe/config.h>
#include <vector>
#include <map>
#include <cstddef>
#include <iosfwd>

typedef int papi_event_t;

class EventSet;

#ifdef HAVE_PAPI

class EventSet
{
private:
  typedef int papi_eventset_t;

  std::vector<papi_event_t> events;
  std::vector<long long> counts;
  papi_eventset_t eventSet;

  void init();
  void handlePAPIError(int error, const char *file, int line) const;

public:
  EventSet();
  EventSet(const papi_event_t *events, size_t count);
  bool isActive();
  void addEvent(papi_event_t event);
  void start();
  void stop();
  void clear();
  void scale(long repetitions);
  std::map<papi_event_t, long long> getCounts() const;
  void printCounts(std::ostream& o, bool verbose) const;
  ~EventSet();
};

#endif

#endif
