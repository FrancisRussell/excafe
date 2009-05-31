#ifndef LIBTRIANGLE_H
#define LIBTRIANGLE_H

// This is copied from triangle.h in libtriangle 1.6 because the Debian modifications expose so much API that
// it breaks the current code. Unfortunately, this might break if the triangle API changes.

extern "C"
{

typedef double TRIREAL;

struct triangulateio {
  TRIREAL *pointlist;                                            /* In / out */
  TRIREAL *pointattributelist;                                   /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  TRIREAL *triangleattributelist;                                /* In / out */
  TRIREAL *trianglearealist;                                      /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  TRIREAL *holelist;                     /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  TRIREAL *regionlist;                   /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  TRIREAL *normlist;             /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);
void trifree(void *memptr);

}
#endif
