//
// Created by Dan Wilkins on 2/5/20.
//

#ifndef RAYTRACE_CPU_HEALPIX_H
#define RAYTRACE_CPU_HEALPIX_H

#define isqrt(x) ((int)sqrt(x+0.5))

// returns the normalized Cartesian coordinate given an x,y coordinate for the given face
template <typename T>
void xyf2loc(T x, T y, int face, T* vec) {
    const static int jrll[] = {2,2,2,2,3,3,3,3,4,4,4,4};
    const static int jpll[] = {1,3,5,7,0,2,4,6,1,3,5,7};

    T z, phi, sinTheta;
    T jr = jrll[face] - x - y;
    T nr;
    if (jr<1) {
        nr = jr;
        T tmp = nr*nr/3.;
        z = 1 - tmp;
    } else if (jr>3) {
        nr = 4-jr;
        T tmp = nr*nr/3.;
        z = tmp - 1;
    } else {
        nr = 1;
        z = (2-jr)*2./3.;
    }
    T tmp=jpll[face]*nr+x-y;
    if (tmp<0) tmp+=8;
    if (tmp>=8) tmp-=8;
    phi = (nr<1e-15) ? 0 : (0.25*M_PI*tmp)/nr;

    // convert to Cartesian
    // NB apply a slight rotation in phi to stop things lining up along axes
    sinTheta = sqrt((1.0-z)*(1.0+z));
    vec[0] = sinTheta*cos(phi+0.05);
    vec[1] = sinTheta*sin(phi+0.05);
    vec[2] = z;
}

// gets the face and pixel indices from one NEST index
template <typename T>
void ring2xyf(int order, int pix, int* ix, int* iy, int* face_num)
{
    const static int jrll[] = {2,2,2,2,3,3,3,3,4,4,4,4};
    const static int jpll[] = {1,3,5,7,0,2,4,6,1,3,5,7};

    int iring, iphi, kshift, nr;

    // back out all the Healpix parameters
    // int order  = grbp->idim[0]; // gives the order of the healpix
    int nside_  = 1<<order;
    int nl2 = 2*nside_;
    int npface_ = nside_<<order;
    int ncap_   = (npface_-nside_)<<1;
    int npix_   = 12*npface_;

    if (pix<ncap_) // North Polar cap
    {
        iring = (1+isqrt(1+2*pix))>>1; //counted from North pole
        iphi  = (pix+1) - 2*iring*(iring-1);
        kshift = 0;
        nr = iring;
        *face_num=(iphi-1)/nr;
    }
    else if (pix<(npix_-ncap_)) // Equatorial region
    {
        int ip = pix - ncap_;
        int tmp = (order>=0) ? ip>>(order+2) : ip/(4*nside_);
        iring = tmp+nside_;
        iphi = ip-tmp*4*nside_ + 1;
        kshift = (iring+nside_)&1;
        nr = nside_;
        int ire = iring-nside_+1,
                irm = nl2+2-ire;
        int ifm = iphi - ire/2 + nside_ -1,
                ifp = iphi - irm/2 + nside_ -1;
        if (order>=0)
        { ifm >>= order; ifp >>= order; }
        else
        { ifm /= nside_; ifp /= nside_; }
        *face_num = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
    else // South Polar cap
    {
        int ip = npix_ - pix;
        iring = (1+isqrt(2*ip-1))>>1; //counted from South pole
        iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
        kshift = 0;
        nr = iring;
        iring = 2*nl2-iring;
        *face_num = 8 + (iphi-1)/nr;
    }

    int irt = iring - (jrll[*face_num]*nside_) + 0;
    int ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
    if (ipt>=nl2) ipt-=8*nside_;

    *ix =  (ipt-irt) >>1;
    *iy = (-ipt-irt) >>1;
}

template <typename T>
void get_pixel_vectors(int order, int pix, T* corners, T* center) {

    // int ii, jj,
    int ix, iy, i, j, k, m, face;
    int nside_, step;
    // T dx[2],  x0[2]
    T xc, yc, dc, d;


    nside_  = 1<<order;
    step  = 1;
    ring2xyf<T>(order, pix, &ix, &iy, &face); // only RING ordering for now
    dc = 0.5 / nside_;
    xc = (ix + 0.5)/nside_, yc = (iy + 0.5)/nside_;
    d = 1.0/(step*nside_);
    for(i=0, j=step, k=2*step, m=3*step; i < step; ++i, ++j, ++k, ++m) {
        xyf2loc<T>(xc+dc-i*d, yc+dc, face, &corners[3*i]);
        xyf2loc<T>(xc-dc, yc+dc-i*d, face, &corners[3*j]);
        xyf2loc<T>(xc-dc+i*d, yc-dc, face, &corners[3*k]);
        xyf2loc<T>(xc+dc, yc-dc+i*d, face, &corners[3*m]);
    }


    // average the center vertex for convenience
    center[0] = 0.25*(corners[3*0+0]+corners[3*1+0]+corners[3*2+0]+corners[3*3+0]);
    center[1] = 0.25*(corners[3*0+1]+corners[3*1+1]+corners[3*2+1]+corners[3*3+1]);
    center[2] = 0.25*(corners[3*0+2]+corners[3*1+2]+corners[3*2+2]+corners[3*3+2]);
}

#endif //RAYTRACE_CPU_HEALPIX_H
