#ifndef EXCAFE_UTIL_ALIGNED_ALLOC_HPP
#define EXCAFE_UTIL_ALIGNED_ALLOC_HPP

#include <cstdlib>

namespace excafe
{

namespace util
{

void *aligned_alloc(size_t alignment, size_t size);

void aligned_free(void *data);

}

}

#endif
