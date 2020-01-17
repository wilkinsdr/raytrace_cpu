#include <iostream>
#include <cmath>
using namespace std;

#include "include/text_output.h"

int main()
{
    const double rsph = 10;
    const double rmin = 0.1;
    const double rstar = 1;
    const double V = 0.2;
    const int Nx = 1000;
    const double dz = 0.01;
    const double en0 = 0.1;
    const double enmax = 1.5;
    const int Nen = 1000;
    const bool logbin_en = false;
    const double dens = 1;
    const double tau = 1.5;
    const double line_emis = 0.000001;

    const double dx = 2 * rsph / Nx;
    const int num_rays = Nx * Nx;

    double *energy, **emis, **absorb, *obs_emis, *obs_continuum, *obs_total;

    const double den = (logbin_en) ? exp(log(enmax / en0) / (Nen - 1)) : (enmax - en0) / (Nen - 1);
    energy = new double[Nen];
    for (int ien = 0; ien < Nen; ien++)
        energy[ien] = (logbin_en) ? en0 * pow(den, ien) : en0 + den * ien;

    emis = new double *[num_rays];
    absorb = new double *[num_rays];

    for (int i = 0; i < num_rays; i++)
    {
        emis[i] = new double[Nen];
        absorb[i] = new double[Nen];
        for (int j = 0; j < Nen; j++)
        {
            emis[i][j] = 0;
            absorb[i][j] = 0;
        }
    }

    for (int ix = 0; ix < Nx; ix++)
    {
        const double x = -1 * rsph + ix * dx;
        for (int iy = 0; iy < Nx; iy++)
        {
            const int ray = ix * Nx + iy;
            const double y = -1 * rsph + iy * dx;

            for (double z = rsph; z > (-1 * rsph); z -= dz)
            {
                const double r = sqrt(x * x + y * y + z * z);
                double this_V = V * sqrt(r) / sqrt(rsph);
                const double los_v = this_V * z / sqrt(x * x + y * y + z * z);
                const double energy = 1 + los_v;
                const int ien = (logbin_en) ? static_cast<int>( log(energy / en0) / log(den)) : static_cast<int>(
                        (energy - en0) / den);

                if (ien >= 0 && ien < Nen)
                {
                    if (r < rsph && r > rmin)
                    {
                        emis[ray][ien] += (1. / (r * r)) * dz * dens * pow(energy, 3);
                        if (z > 0) absorb[ray][ien] += dz * dens;
                    }
                }
            }
        }
    }

    const int ray0 = (rsph / dx) * Nx + (rsph / dx);

    obs_emis = new double[Nen];

    double emis_sum = 0;

    for (int ien = 0; ien < Nen; ien++)
    {
        obs_emis[ien] = 0;
        for (int ray = 0; ray < num_rays; ray++)
            obs_emis[ien] += emis[ray][ien];
        emis_sum += obs_emis[ien];
    }

    obs_continuum = new double[Nen];

    double tau_total = 0;

    for (int ien = 0; ien < Nen; ien++)
    {
        tau_total += absorb[ray0][ien];
    }

    double continuum_sum = 0;

//    for(int ien=0; ien<Nen; ien++)
//    {
//        obs_continuum[ien] = 1;
//        obs_continuum[ien] *= exp(-1*(tau/tau_total)*absorb[ray0][ien]);
//        continuum_sum += obs_continuum[ien];
//    }

    for(int ien=0; ien<Nen; ien++)
        obs_continuum[ien] = 0;

    const int istart = static_cast<int>((-rstar + rsph) / dx);
    const int iend = static_cast<int>((rstar + rsph) / dx);
    for(int ix=istart; ix<=iend; ix++)
    {
        const double x = -1*rsph + ix*dx;
        for(int iy=istart; iy<=iend; iy++)
        {
            const double y = -1*rsph + iy*dx;
            const int ray = ix*Nx + iy;

            if((x*x + y*y) < (rstar*rstar))
                for(int ien=0; ien<Nen; ien++)
                {
                    obs_continuum[ien] += exp(-1*(tau/tau_total)*absorb[ray][ien]);
                    continuum_sum += obs_continuum[ien];
                }
        }
    }

    obs_total = new double[Nen];
    for(int ien=0; ien<Nen; ien++)
    {
        obs_total[ien] = (line_emis/emis_sum) * obs_emis[ien] + (1./continuum_sum) * obs_continuum[ien];
    }

    TextOutput outfile("../dat/pcyg.dat");
    for(int ien=0; ien<Nen; ien++)
    {
        outfile << energy[ien] << obs_emis[ien] << obs_continuum[ien] << obs_total[ien] << endl;
    }


    return 0;
}