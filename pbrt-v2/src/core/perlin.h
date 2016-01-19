#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_PERLIN_H
#define PBRT_CORE_PERLIN_H

// core/splines.h*
#include "pbrt.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <core/splines.h>

class Perlin {
public:
    virtual ~Perlin(){}
    Perlin(float _persistence, float _frequency, int _level, unsigned int seed)
        :persistence(_persistence), frequency(_frequency), level(_level)
    {
        srand(seed);
        for(int i=0;i<level;i++){
            float l_frequency=frequency/pow(2.0f,(float)i);
            float altitude=pow(persistence,(float)i);
            Catmull_Rom spline(true);
            for(int j=0;j<1.0/l_frequency+2;j++){
                spline.Add_Sample((float)j*l_frequency,altitude*(float)rand()/(float)RAND_MAX);
            }
            spline.Build();
            noises.push_back(spline);
        }
    }

    float Evaluate(float x){
        float value=0;
        for(int i=0;i<level;i++){
            value+=noises[i].Evaluate(x);
        }
        return value;
    }
    std::vector<Catmull_Rom> noises;
    float persistence;
    float frequency;
    int level;
};

#endif // PBRT_CORE_SPLINES_H
