#include <excafe/config.h>
#include "event_set.hpp"
#include <papi.h>
#include <map>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <boost/foreach.hpp>
#include <ostream>

#include <iostream>

#define PAPI_CHECK(err) handlePAPIError(err, __FILE__, __LINE__);

void EventSet::handlePAPIError(const int error, const char *const file, const int line) const
{
  if (error < 0)
  {
    std::ostringstream stream;

    stream << "Error received from PAPI in " << file << " at line " << line << ": ";
    const char *const errorString = PAPI_strerror(error);

    if (errorString == NULL)
      stream << "Unrecognised error code (" << error << ")";
    else
      stream << errorString;

    stream << ".";
    throw std::runtime_error(stream.str());
  }
}

void EventSet::init()
{
  if (!PAPI_is_initialized())
  {
    PAPI_CHECK(PAPI_library_init(PAPI_VER_CURRENT));
    PAPI_CHECK(PAPI_multiplex_init());
  }

  eventSet = PAPI_NULL;
  PAPI_CHECK(PAPI_create_eventset(&eventSet));
}

EventSet::EventSet()
{
  init();
}

EventSet::EventSet(const papi_event_t *const _events, const size_t count) :
  events(_events, _events + count), counts(count, 0)
{
  init();
  PAPI_CHECK(PAPI_add_events(eventSet, &events[0], events.size()));
}

bool EventSet::isActive()
{
  int status;
  PAPI_CHECK(PAPI_state(eventSet, &status));

  return (status & PAPI_RUNNING) != 0;
}

void EventSet::addEvent(const papi_event_t event)
{
  PAPI_CHECK(PAPI_add_event(eventSet, event));
  events.push_back(event);
  counts.push_back(0);
}

void EventSet::start()
{
  PAPI_CHECK(PAPI_start(eventSet));
}

void EventSet::stop()
{
  PAPI_CHECK(PAPI_stop(eventSet, &counts[0]));
}

void EventSet::clear()
{
  std::fill(counts.begin(), counts.end(), 0);
}

void EventSet::scale(long repetitions)
{
  BOOST_FOREACH(long long& value, counts)
    value /= repetitions;
}

std::map<papi_event_t, long long> EventSet::getCounts() const
{
  assert(events.size() == counts.size());
  std::map<papi_event_t, long long> countMap;

  for(size_t i = 0; i < events.size(); ++i)
    countMap[events[i]] = counts[i];

  return countMap;
}

void EventSet::printCounts(std::ostream& o, const bool verbose) const
{
  typedef std::pair<papi_event_t, long long> event_freq_t;
  const std::map<papi_event_t, long long> countMap(getCounts());

  BOOST_FOREACH(const event_freq_t& eventFreq, countMap)
  {
    PAPI_event_info_t eventInfo;
    PAPI_CHECK(PAPI_get_event_info(eventFreq.first, &eventInfo));

    o << (verbose ? eventInfo.long_descr : eventInfo.short_descr);
    o << ": " << eventFreq.second << std::endl;
  }
}

EventSet::~EventSet()
{
  if (isActive())
    stop();

  PAPI_CHECK(PAPI_cleanup_eventset(eventSet));
  PAPI_CHECK(PAPI_destroy_eventset(&eventSet));
}
