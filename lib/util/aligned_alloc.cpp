#include <excafe/util/aligned_alloc.hpp>
#include <cstdlib>
#include <stdint.h>

namespace excafe
{

namespace util
{

void *aligned_alloc(const size_t alignment, const size_t size)
{
  char *const data = new char[size + sizeof(size_t) + alignment - 1];
  char *start = data + sizeof(size_t);
  const uintptr_t start_alignment = (uintptr_t) start & (alignment - 1);

  if (start_alignment != 0)
    start += (alignment - start_alignment);

  const size_t offset = start - data;
  ((size_t*) start)[-1] = offset;

  return start;
}

void aligned_free(void *data)
{
  const size_t offset = ((size_t*) data)[-1];
  delete ((char*)data - offset);
}

}

}
