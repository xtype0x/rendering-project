#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_PERLIN2D_H
#define PBRT_CORE_PERLIN2D_H

// core/splines.h*
#include "pbrt.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>

class Perlin2D {
public:
    virtual ~Perlin2D(){}
    Perlin2D(float _persistence, float _frequency, int _level, unsigned int seed)
        :persistence(_persistence), frequency(_frequency), level(_level)
    {
        srand(seed);
        noises.resize(level);
        l_frequencies.resize(level);
        for(int i=0;i<level;i++){
            float l_frequency=frequency/pow(2.0f,(float)i);
            l_frequencies[i]=l_frequency;
            size_t n_samples=1.0/l_frequency+2;
            float altitude=pow(persistence,(float)i);
            noises[i].resize(n_samples);
            for(int j=0;j<n_samples;j++){
                noises[i][j].resize(n_samples);
                for(int k=0;k<n_samples;k++){
                    noises[i][j][k]=altitude*(float)rand()/(float)RAND_MAX;
                }
            }
        }
    }

    float Interpolate(float v1, float v2, float frac) const
    {
        float ft=frac*3.1415926;
        float f=(1-cos(ft))*0.5f;
        return v1*(1.0f-f)+v2*f;
        //return v1*(1.0f-frac)+v2*frac;
    }

    float InterpolateNoise(float x, float y, int i) const
    {
        int index_x=x/l_frequencies[i];
        int index_y=y/l_frequencies[i];
        float frac_x=(x-(float)index_x*l_frequencies[i])/l_frequencies[i];
        float frac_y=(y-(float)index_y*l_frequencies[i])/l_frequencies[i];
        float n1=Interpolate(noises[i][index_x][index_y],noises[i][index_x+1][index_y],frac_x);
        float n2=Interpolate(noises[i][index_x][index_y+1],noises[i][index_x+1][index_y+1],frac_x);
        return Interpolate(n1,n2,frac_y);
    }

    float Evaluate(float x, float y) const
    {
        float value=0;
        for(int i=0;i<level;i++){
            value+=InterpolateNoise(x,y,i);
        }
        return value;
    }
    std::vector<std::vector<std::vector<float> > > noises;
    std::vector<float> l_frequencies;
    float persistence;
    float frequency;
    int level;
};

#endif // PBRT_CORE_SPLINES_H
