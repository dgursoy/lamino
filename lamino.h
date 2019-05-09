#ifndef _lamino_h
#define _lamino_h

#include "string.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327
#endif


void project(
    float* obj, 
    const int ox,  
    const int oy,  
    const int oz, 
    const float* gridx, 
    const float* gridy, 
    const float* gridz,
    const float* phi,
    const float* theta,
    const float* detgridx,
    const float* detgridy,
    const int dx,  
    const int dy,  
    float* proj, 
    const int px,  
    const int py,  
    const int pz);


void art(
    float* prj, 
    const int px, 
    const int py, 
    const int pz,
    const float* phi,
    const float* theta,
    const float* detgridx,
    const float* detgridy,
    const int dx,
    const int dy,
    float* rec,
    const int ox, 
    const int oy, 
    const int oz,
    const float* gridx, 
    const float* gridy, 
    const float* gridz,
    const int iters);


void mlem(
    float* prj, 
    const int px, 
    const int py, 
    const int pz,
    const float* phi,
    const float* theta,
    const float* detgridx,
    const float* detgridy,
    const int dx,
    const int dy,
    float* rec,
    const int ox, 
    const int oy, 
    const int oz,
    const float* gridx, 
    const float* gridy, 
    const float* gridz,
    const int iters);

#endif
