#include "lamino.h"
#include <time.h>

void
project(
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
    const int pz)
{
    int m, n, mm, nn, tt, a, indproj;
    float srcx, srcy, srcz, detx, dety, detz;
    float _srcx, _srcy, _srcz, _detx, _dety, _detz;
    float t, temp1, temp2, temp3;
    int numcoord, numind;
    float _dx, _dy, _dz;
    int gx = ox + 1;
    int gy = oy + 1;
    int gz = oz + 1;
    float* tlen = calloc(gx+gy+gz, sizeof(float));
    float* coordx = calloc(gx+gy+gz, sizeof(float));
    float* coordy = calloc(gx+gy+gz, sizeof(float));
    float* coordz = calloc(gx+gy+gz, sizeof(float));
    float* dist = calloc(gx+gy+gz, sizeof(float));
    float* mx = calloc(gx+gy+gz, sizeof(float));
    float* my = calloc(gx+gy+gz, sizeof(float));
    float* mz = calloc(gx+gy+gz, sizeof(float));
    int* ind = calloc(gx+gy+gz, sizeof(int));
    int* ix = calloc(gx+gy+gz, sizeof(int));
    int* iy = calloc(gx+gy+gz, sizeof(int));
    int* iz = calloc(gx+gy+gz, sizeof(int));

    clock_t start, end;
    double cpu_time_used1 = 0;
    double cpu_time_used2 = 0; 
    double cpu_time_used3 = 0; 
    double cpu_time_used4 = 0; 
    double cpu_time_used5 = 0; 
    double cpu_time_used6 = 0; 
    double cpu_time_used7 = 0; 

    _srcz = -0.5;
    _detz = 0.5;

    for (tt = 0; tt < px; tt++)
    {
        printf("%d\n", tt);

        // Make sure the angles are in [0, 360) degrees
        float rot = fmodf(theta[tt], 2.0f * (float) M_PI);
        float tilt = fmodf(phi[tt], 2.0f * (float) M_PI);

        // Cache trigonometric coefficients
        float cost, sint, cosr, sinr, sinsin, sincos, cossin, coscos;
        cost = cosf(tilt);
        sint = sinf(tilt);
        cosr = cosf(rot);
        sinr = sinf(rot);
        sinsin = sint * sinr;
        sincos = sint * cosr;
        cossin = cost * sinr;
        coscos = cost * cosr;

        for (mm = 0; mm < dy; mm++)
        {
            _dety = detgridy[mm];
            _srcy = _dety;

            for (nn = 0; nn < dx; nn++)
            {
                _detx = detgridx[nn];
                _srcx = _detx;

                // Tilt and rotate
                srcx = _srcx * cost - _srcz * sint;
                srcy = _srcy * cosr - _srcx * sinsin - _srcz * cossin;
                srcz = _srcy * sinr + _srcx * sincos + _srcz * coscos;
                detx = _detx * cost - _detz * sint;
                dety = _dety * cosr - _detx * sinsin - _detz * cossin;
                detz = _dety * sinr + _detx * sincos + _detz * coscos;

                // Calculate plane intersections
                start = clock();
                numcoord = 0;
                for(m = 0; m < gx; m++)
                { 
                    t = (gridx[m] - srcx) / (detx - srcx);
                    temp1 = srcy + t * (dety - srcy);
                    temp2 = srcz + t * (detz - srcz);
                    if (temp1 >= gridy[0])
                    {
                        if  (temp1 <= gridy[oy])
                        {
                            if (temp2 >= gridz[0])
                            {
                                if  (temp2 <= gridz[oz])
                                {
                                    tlen[numcoord] = t;
                                    coordx[numcoord] = gridx[m];
                                    coordy[numcoord] = temp1;
                                    coordz[numcoord] = temp2;
                                    numcoord++;
                                }
                            }
                        }
                    }
                }
                for(m = 0; m < gy; m++)
                {
                    t = (gridy[m] - srcy) / (dety - srcy);
                    temp1 = srcx + t * (detx - srcx);
                    temp2 = srcz + t * (detz - srcz);
                    if (temp1 >= gridx[0])
                    {
                        if  (temp1 <= gridx[ox])
                        {
                            if (temp2 >= gridz[0])
                            {
                                if  (temp2 <= gridz[oz])
                                {
                                    tlen[numcoord] = t;
                                    coordx[numcoord] = temp1;
                                    coordy[numcoord] = gridy[m];
                                    coordz[numcoord] = temp2;
                                    numcoord++;
                                }
                            }
                        }
                    }
                }
                for(m = 0; m < gz; m++)
                {
                    t = (gridz[m] - srcz) / (detz - srcz);
                    temp1 = srcx + t * (detx - srcx);
                    temp2 = srcy + t * (dety - srcy);
                    if (temp1 >= gridx[0])
                    {
                        if  (temp1 <= gridx[ox])
                        {
                            if (temp2 >= gridy[0])
                            {
                                if  (temp2 <= gridy[oy])
                                {
                                    tlen[numcoord] = t;
                                    coordx[numcoord] = temp1;
                                    coordy[numcoord] = temp2;
                                    coordz[numcoord] = gridz[m];
                                    numcoord++;
                                }
                            }
                        }
                    }
                }
                end = clock();
                cpu_time_used1 += ((double) (end - start)) / CLOCKS_PER_SEC;

                // Sort coordinates in ascending order according to tlen
                start = clock();
                for (m = 0; m < numcoord; m++)
                {
                    for (n = 0; n < numcoord; n++)
                    {
                        if (tlen[n] > tlen[m])
                        {
                            temp1 = tlen[m];
                            tlen[m] = tlen[n];
                            tlen[n] = temp1;

                            temp2 = coordx[m];
                            coordx[m] = coordx[n];
                            coordx[n] = temp2;

                            temp3 = coordy[m];
                            coordy[m] = coordy[n];
                            coordy[n] = temp3;

                            temp3 = coordz[m];
                            coordz[m] = coordz[n];
                            coordz[n] = temp3;
                        }  
                    }
                }   
                end = clock();
                cpu_time_used2 += ((double) (end - start)) / CLOCKS_PER_SEC;

                // Calculate ray path lengths in each voxel
                start = clock();
                numind = numcoord-1;
                for (m = 0; m < numind; m++)
                {
                    _dx = (coordx[m + 1] - coordx[m]) * (coordx[m + 1] - coordx[m]);
                    _dy = (coordy[m + 1] - coordy[m]) * (coordy[m + 1] - coordy[m]);
                    _dz = (coordz[m + 1] - coordz[m]) * (coordz[m + 1] - coordz[m]);
                    dist[m] = sqrtf(_dx + _dy + _dz);
                }
                end = clock();
                cpu_time_used3 += ((double) (end - start)) / CLOCKS_PER_SEC;

                // Calculate middle points of ray segments in each voxel
                start = clock();
                for (m = 0; m < numind; m++)
                {
                    mx[m] = coordx[m] + 0.5 * (coordx[m + 1] - coordx[m]);
                    my[m] = coordy[m] + 0.5 * (coordy[m + 1] - coordy[m]);
                    mz[m] = coordz[m] + 0.5 * (coordz[m + 1] - coordz[m]);
                }
                end = clock();
                cpu_time_used4 += ((double) (end - start)) / CLOCKS_PER_SEC;

                // Calc object's indices of ray path segments
                start = clock();
                for (m = 0; m < numind; m++)
                {
                    a = 0;
                    for (n = 1; n < gx; n++)
                    {
                        if (mx[m] <= gridx[n])
                        {
                            ix[m] = a;
                            break;
                        }
                        a++;
                    }
                }
                for (m = 0; m < numind; m++)
                {
                    a = 0;
                    for (n = 1; n < gy; n++)
                    {
                        if (my[m] <= gridy[n])
                        {
                            iy[m] = a;
                            break;
                        }
                        a++;
                    }
                }
                for (m = 0; m < numind; m++)
                {
                    a = 0;
                    for (n = 1; n < gz; n++)
                    {
                        if (mz[m] <= gridz[n])
                        {
                            iz[m] = a;
                            break;
                        }
                        a++;
                    }
                }
                end = clock();
                cpu_time_used5 += ((double) (end - start)) / CLOCKS_PER_SEC;

                start = clock();
                for (m = 0; m < numind; m++)
                {
                    ind[m] = iz[m] + oz * iy[m] + oy * oz * ix[m];
                }
                end = clock();
                cpu_time_used6 += ((double) (end - start)) / CLOCKS_PER_SEC;
                
                // Calculate ray sum
                start = clock();
                indproj = mm + pz * nn + py * pz * tt;
                for (m = 0; m < numind; m++)
                {
                    proj[indproj] += obj[ind[m]] * dist[m];
                }
                end = clock();
                cpu_time_used7 += ((double) (end - start)) / CLOCKS_PER_SEC;
            }
        }
    }
                
    printf ("cpu_time_used1=%f \n", cpu_time_used1);
    printf ("cpu_time_used2=%f \n", cpu_time_used2);
    printf ("cpu_time_used3=%f \n", cpu_time_used3);
    printf ("cpu_time_used4=%f \n", cpu_time_used4);
    printf ("cpu_time_used5=%f \n", cpu_time_used5);
    printf ("cpu_time_used6=%f \n", cpu_time_used6);
    printf ("cpu_time_used7=%f \n", cpu_time_used7);

    free(tlen);
    free(coordx);
    free(coordy);
    free(coordz);
    free(dist);
    free(mx);
    free(my);
    free(mz);
    free(ind);
    free(ix);
    free(iy);
    free(iz);
}


void
art(
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
    const int iters)
{
    int m, n, ii, mm, nn, tt, a, indproj;
    float srcx, srcy, srcz, detx, dety, detz;
    float _srcx, _srcy, _srcz, _detx, _dety, _detz;
    float t, temp1, temp2, temp3;
    int numcoord, numind;
    float _dx, _dy, _dz;
    int gx = ox + 1;
    int gy = oy + 1;
    int gz = oz + 1;
    float sim;
    float* tlen = calloc(gx+gy+gz, sizeof(float));
    float* coordx = calloc(gx+gy+gz, sizeof(float));
    float* coordy = calloc(gx+gy+gz, sizeof(float));
    float* coordz = calloc(gx+gy+gz, sizeof(float));
    float* dist = calloc(gx+gy+gz, sizeof(float));
    float* mx = calloc(gx+gy+gz, sizeof(float));
    float* my = calloc(gx+gy+gz, sizeof(float));
    float* mz = calloc(gx+gy+gz, sizeof(float));
    int* ind = calloc(gx+gy+gz, sizeof(int));
    int* ix = calloc(gx+gy+gz, sizeof(int));
    int* iy = calloc(gx+gy+gz, sizeof(int));
    int* iz = calloc(gx+gy+gz, sizeof(int));

    _srcz = -0.5;
    _detz = 0.5;

    for (ii = 0; ii < iters; ii++)
    {
        for (tt = 0; tt < px; tt++)
        {
            printf("%d\n", tt);

            // Make sure the angles are in [0, 360) degrees
            float rot = fmodf(theta[tt], 2.0f * (float) M_PI);
            float tilt = fmodf(phi[tt], 2.0f * (float) M_PI);

            // Cache trigonometric coefficients
            float cost, sint, cosr, sinr, sinsin, sincos, cossin, coscos;
            cost = cosf(tilt);
            sint = sinf(tilt);
            cosr = cosf(rot);
            sinr = sinf(rot);
            sinsin = sint * sinr;
            sincos = sint * cosr;
            cossin = cost * sinr;
            coscos = cost * cosr;

            for (mm = 0; mm < dy; mm++)
            {
                _dety = detgridy[mm];
                _srcy = _dety;

                for (nn = 0; nn < dx; nn++)
                {
                    _detx = detgridx[nn];
                    _srcx = _detx;

                    // Tilt and rotate
                    srcx = _srcx * cost - _srcz * sint;
                    srcy = _srcy * cosr - _srcx * sinsin - _srcz * cossin;
                    srcz = _srcy * sinr + _srcx * sincos + _srcz * coscos;
                    detx = _detx * cost - _detz * sint;
                    dety = _dety * cosr - _detx * sinsin - _detz * cossin;
                    detz = _dety * sinr + _detx * sincos + _detz * coscos;

                    // Calculate plane intersections
                    numcoord = 0;
                    for(m = 0; m < gx; m++)
                    { 
                        t = (gridx[m] - srcx) / (detx - srcx);
                        temp1 = srcy + t * (dety - srcy);
                        temp2 = srcz + t * (detz - srcz);
                        if (temp1 >= gridy[0])
                        {
                            if  (temp1 <= gridy[oy])
                            {
                                if (temp2 >= gridz[0])
                                {
                                    if  (temp2 <= gridz[oz])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = gridx[m];
                                        coordy[numcoord] = temp1;
                                        coordz[numcoord] = temp2;
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }
                    for(m = 0; m < gy; m++)
                    {
                        t = (gridy[m] - srcy) / (dety - srcy);
                        temp1 = srcx + t * (detx - srcx);
                        temp2 = srcz + t * (detz - srcz);
                        if (temp1 >= gridx[0])
                        {
                            if  (temp1 <= gridx[ox])
                            {
                                if (temp2 >= gridz[0])
                                {
                                    if  (temp2 <= gridz[oz])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = temp1;
                                        coordy[numcoord] = gridy[m];
                                        coordz[numcoord] = temp2;
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }
                    for(m = 0; m < gz; m++)
                    {
                        t = (gridz[m] - srcz) / (detz - srcz);
                        temp1 = srcx + t * (detx - srcx);
                        temp2 = srcy + t * (dety - srcy);
                        if (temp1 >= gridx[0])
                        {
                            if  (temp1 <= gridx[ox])
                            {
                                if (temp2 >= gridy[0])
                                {
                                    if  (temp2 <= gridy[oy])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = temp1;
                                        coordy[numcoord] = temp2;
                                        coordz[numcoord] = gridz[m];
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }

                    // Sort coordinates in ascending order according to tlen
                    for (m = 0; m < numcoord; m++)
                    {
                        for (n = 0; n < numcoord; n++)
                        {
                            if (tlen[n] > tlen[m])
                            {
                                temp1 = tlen[m];
                                tlen[m] = tlen[n];
                                tlen[n] = temp1;

                                temp2 = coordx[m];
                                coordx[m] = coordx[n];
                                coordx[n] = temp2;

                                temp3 = coordy[m];
                                coordy[m] = coordy[n];
                                coordy[n] = temp3;

                                temp3 = coordz[m];
                                coordz[m] = coordz[n];
                                coordz[n] = temp3;
                            }  
                        }
                    }

                    // Calculate ray path lengths in each voxel
                    numind = numcoord-1;
                    for (m = 0; m < numind; m++)
                    {
                        _dx = (coordx[m + 1] - coordx[m]) * (coordx[m + 1] - coordx[m]);
                        _dy = (coordy[m + 1] - coordy[m]) * (coordy[m + 1] - coordy[m]);
                        _dz = (coordz[m + 1] - coordz[m]) * (coordz[m + 1] - coordz[m]);
                        dist[m] = sqrtf(_dx + _dy + _dz);
                    }

                    // Calculate middle points of ray segments in each voxel
                    for (m = 0; m < numind; m++)
                    {
                        mx[m] = coordx[m] + 0.5 * (coordx[m + 1] - coordx[m]);
                        my[m] = coordy[m] + 0.5 * (coordy[m + 1] - coordy[m]);
                        mz[m] = coordz[m] + 0.5 * (coordz[m + 1] - coordz[m]);
                    }

                    // Calc object's indices of ray path segments
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gx; n++)
                        {
                            if (mx[m] <= gridx[n])
                            {
                                ix[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gy; n++)
                        {
                            if (my[m] <= gridy[n])
                            {
                                iy[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gz; n++)
                        {
                            if (mz[m] <= gridz[n])
                            {
                                iz[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        ind[m] = iz[m] + oz * iy[m] + oy * oz * ix[m];
                    }
                    
                    // Calculate ray sum
                    sim = 0;
                    for (m = 0; m < numind; m++)
                    {
                        sim += rec[ind[m]] * dist[m];
                    }

                    // Update
                    float sum_dist2 = 0;
                    for (m = 0; m < numind; m++)
                    {
                        sum_dist2 += dist[m] * dist[m];
                    }
                    if(sum_dist2 != 0.0f)
                    {
                        indproj = mm + pz * nn + py * pz * tt;
                        float upd = (prj[indproj] - sim) / sum_dist2;
                        for(m = 0; m < numind; m++)
                        {
                            rec[ind[m]] += upd * dist[m];
                        }
                    }
                }
            }
        }
    }

    free(tlen);
    free(coordx);
    free(coordy);
    free(coordz);
    free(dist);
    free(mx);
    free(my);
    free(mz);
    free(ind);
    free(ix);
    free(iy);
    free(iz);
}


void
mlem(
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
    const int iters)
{
    int m, n, ii, mm, nn, tt, a, indproj;
    float srcx, srcy, srcz, detx, dety, detz;
    float _srcx, _srcy, _srcz, _detx, _dety, _detz;
    float t, temp1, temp2, temp3;
    int numcoord, numind;
    float _dx, _dy, _dz;
    int gx = ox + 1;
    int gy = oy + 1;
    int gz = oz + 1;
    float sim;
    float* tlen = calloc(gx+gy+gz, sizeof(float));
    float* coordx = calloc(gx+gy+gz, sizeof(float));
    float* coordy = calloc(gx+gy+gz, sizeof(float));
    float* coordz = calloc(gx+gy+gz, sizeof(float));
    float* dist = calloc(gx+gy+gz, sizeof(float));
    float* mx = calloc(gx+gy+gz, sizeof(float));
    float* my = calloc(gx+gy+gz, sizeof(float));
    float* mz = calloc(gx+gy+gz, sizeof(float));
    int* ind = calloc(gx+gy+gz, sizeof(int));
    int* ix = calloc(gx+gy+gz, sizeof(int));
    int* iy = calloc(gx+gy+gz, sizeof(int));
    int* iz = calloc(gx+gy+gz, sizeof(int));
    float sumdist2, upd;
    float* sumdist = calloc(ox+oy+oz, sizeof(float));
    float* update =  calloc(ox+oy+oz, sizeof(float));

    _srcz = -0.5;
    _detz = 0.5;

    for (ii = 0; ii < iters; ii++)
    {
        for (tt = 0; tt < px; tt++)
        {
            printf("%d\n", tt);

            // Make sure the angles are in [0, 360) degrees
            float rot = fmodf(theta[tt], 2.0f * (float) M_PI);
            float tilt = fmodf(phi[tt], 2.0f * (float) M_PI);

            // Cache trigonometric coefficients
            float cost, sint, cosr, sinr, sinsin, sincos, cossin, coscos;
            cost = cosf(tilt);
            sint = sinf(tilt);
            cosr = cosf(rot);
            sinr = sinf(rot);
            sinsin = sint * sinr;
            sincos = sint * cosr;
            cossin = cost * sinr;
            coscos = cost * cosr;

            for (mm = 0; mm < dy; mm++)
            {
                _dety = detgridy[mm];
                _srcy = _dety;

                for (nn = 0; nn < dx; nn++)
                {
                    _detx = detgridx[nn];
                    _srcx = _detx;

                    // Tilt and rotate
                    srcx = _srcx * cost - _srcz * sint;
                    srcy = _srcy * cosr - _srcx * sinsin - _srcz * cossin;
                    srcz = _srcy * sinr + _srcx * sincos + _srcz * coscos;
                    detx = _detx * cost - _detz * sint;
                    dety = _dety * cosr - _detx * sinsin - _detz * cossin;
                    detz = _dety * sinr + _detx * sincos + _detz * coscos;

                    // Calculate plane intersections
                    numcoord = 0;
                    for(m = 0; m < gx; m++)
                    { 
                        t = (gridx[m] - srcx) / (detx - srcx);
                        temp1 = srcy + t * (dety - srcy);
                        temp2 = srcz + t * (detz - srcz);
                        if (temp1 >= gridy[0])
                        {
                            if  (temp1 <= gridy[oy])
                            {
                                if (temp2 >= gridz[0])
                                {
                                    if  (temp2 <= gridz[oz])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = gridx[m];
                                        coordy[numcoord] = temp1;
                                        coordz[numcoord] = temp2;
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }
                    for(m = 0; m < gy; m++)
                    {
                        t = (gridy[m] - srcy) / (dety - srcy);
                        temp1 = srcx + t * (detx - srcx);
                        temp2 = srcz + t * (detz - srcz);
                        if (temp1 >= gridx[0])
                        {
                            if  (temp1 <= gridx[ox])
                            {
                                if (temp2 >= gridz[0])
                                {
                                    if  (temp2 <= gridz[oz])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = temp1;
                                        coordy[numcoord] = gridy[m];
                                        coordz[numcoord] = temp2;
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }
                    for(m = 0; m < gz; m++)
                    {
                        t = (gridz[m] - srcz) / (detz - srcz);
                        temp1 = srcx + t * (detx - srcx);
                        temp2 = srcy + t * (dety - srcy);
                        if (temp1 >= gridx[0])
                        {
                            if  (temp1 <= gridx[ox])
                            {
                                if (temp2 >= gridy[0])
                                {
                                    if  (temp2 <= gridy[oy])
                                    {
                                        tlen[numcoord] = t;
                                        coordx[numcoord] = temp1;
                                        coordy[numcoord] = temp2;
                                        coordz[numcoord] = gridz[m];
                                        numcoord++;
                                    }
                                }
                            }
                        }
                    }

                    // Sort coordinates in ascending order according to tlen
                    for (m = 0; m < numcoord; m++)
                    {
                        for (n = 0; n < numcoord; n++)
                        {
                            if (tlen[n] > tlen[m])
                            {
                                temp1 = tlen[m];
                                tlen[m] = tlen[n];
                                tlen[n] = temp1;

                                temp2 = coordx[m];
                                coordx[m] = coordx[n];
                                coordx[n] = temp2;

                                temp3 = coordy[m];
                                coordy[m] = coordy[n];
                                coordy[n] = temp3;

                                temp3 = coordz[m];
                                coordz[m] = coordz[n];
                                coordz[n] = temp3;
                            }  
                        }
                    }

                    // Calculate ray path lengths in each voxel
                    numind = numcoord-1;
                    for (m = 0; m < numind; m++)
                    {
                        _dx = (coordx[m + 1] - coordx[m]) * (coordx[m + 1] - coordx[m]);
                        _dy = (coordy[m + 1] - coordy[m]) * (coordy[m + 1] - coordy[m]);
                        _dz = (coordz[m + 1] - coordz[m]) * (coordz[m + 1] - coordz[m]);
                        dist[m] = sqrtf(_dx + _dy + _dz);
                    }

                    // Calculate middle points of ray segments in each voxel
                    for (m = 0; m < numind; m++)
                    {
                        mx[m] = coordx[m] + 0.5 * (coordx[m + 1] - coordx[m]);
                        my[m] = coordy[m] + 0.5 * (coordy[m + 1] - coordy[m]);
                        mz[m] = coordz[m] + 0.5 * (coordz[m + 1] - coordz[m]);
                    }

                    // Calc object's indices of ray path segments
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gx; n++)
                        {
                            if (mx[m] <= gridx[n])
                            {
                                ix[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gy; n++)
                        {
                            if (my[m] <= gridy[n])
                            {
                                iy[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        a = 0;
                        for (n = 1; n < gz; n++)
                        {
                            if (mz[m] <= gridz[n])
                            {
                                iz[m] = a;
                                break;
                            }
                            a++;
                        }
                    }
                    for (m = 0; m < numind; m++)
                    {
                        ind[m] = iz[m] + oz * iy[m] + oy * oz * ix[m];
                    }
                    
                    // Calculate ray sum
                    sim = 0;
                    for (m = 0; m < numind; m++)
                    {
                        sim += rec[ind[m]] * dist[m];
                    }

                    // Update
                    sumdist2 = 0;
                    for (m = 0; m < numind; m++)
                    {
                        sumdist[ind[m]] += dist[m];
                        sumdist2 += dist[m] * dist[m];
                    }
                    if (sumdist2 != 0.0f)
                    {
                        indproj = mm + pz * nn + py * pz * tt;
                        upd = prj[indproj] / sim;
                        for (m = 0; m < numind; m++)
                        {
                            update[ind[m]] += upd * dist[m];
                        }
                    }
                }
            }
        }

        for(m = 0; m < ox * oy * oz; m++)
        {
            if(sumdist[m] != 0.0f)
            {
                rec[m] *= update[m] / sumdist[m];
            }
        }
    }

    free(tlen);
    free(coordx);
    free(coordy);
    free(coordz);
    free(dist);
    free(mx);
    free(my);
    free(mz);
    free(ind);
    free(ix);
    free(iy);
    free(iz);
    free(sumdist);
    free(update);
}


