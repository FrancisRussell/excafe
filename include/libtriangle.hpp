#ifndef LIBTRIANGLE_HPP
#define LIBTRIANGLE_HPP

// This is copied from triangle.h in libtriangle 1.6 because the Debian modifications expose so much API that
// it breaks the current code. Unfortunately, this might break if the triangle API changes.

extern "C"
{

// We assume REAL was defined as double

struct triangulateio {
  double *pointlist;                                             /* In / out */
  double *pointattributelist;                                    /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  double *triangleattributelist;                                 /* In / out */
  double *trianglearealist;                                       /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  double *holelist;                      /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  double *regionlist;                    /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  double *normlist;              /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

void triangulate(char *triswitches, struct triangulateio *in, struct triangulateio *out, struct triangulateio *vorout);
void trifree(void *memptr);

}
#endif
