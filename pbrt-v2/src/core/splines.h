#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_SPLINES_H
#define PBRT_CORE_SPLINES_H

// core/splines.h*
#include "pbrt.h"
#include <vector>
#include <iostream>
#include <algorithm>

class Catmull_Rom {
public:
    virtual ~Catmull_Rom(){}
    Catmull_Rom(bool _even_sorted_sample=false)
        :even_sorted_sample(_even_sorted_sample)
    {}
    void Add_Sample(float x, float y){
        samples.push_back(std::pair<float,float>(x,y));
    }

    static bool SortSample (const std::pair<float,float>& i,const std::pair<float,float>& j) { return (i.first<j.first); }
    void Build(){
        if(samples.size()<2) return;
        if(!even_sorted_sample)
            std::sort(samples.begin(),samples.end(),SortSample);
        else
            step_size=(samples[samples.size()-1].first-samples[0].first)/(float)(samples.size()-1);
        tangents.resize(samples.size());
        tangents[0]=(samples[1].second-samples[0].second)/(samples[1].first-samples[0].first);
        tangents[tangents.size()-1]=(samples[tangents.size()-1].second-samples[tangents.size()-2].second)/(samples[tangents.size()-1].first-samples[tangents.size()-2].first);
        for(size_t i=1;i<tangents.size()-1;i++)
            tangents[i]=(samples[i+1].second-samples[i-1].second)/(samples[i+1].first-samples[i-1].first);
    }

    float Evaluate(float x)
    {
        if(x<samples[0].first) return samples[0].second;
        if(x>samples[samples.size()-1].first) return samples[samples.size()-1].second;
        size_t first=0, last=samples.size()-1, midPoint;
        if(!even_sorted_sample){
            while(first+1<last){
                midPoint = (first+last)/2;
                if (x < samples[midPoint].first)
                    last = midPoint;
                else {
                    first = midPoint;
                    break;
                }
            }
            assert(x>=samples[first].first&&x<samples[first+1].first);
        }
        else{
            first=(x-samples[first].first)/step_size;
        }
        float t=(x-samples[first].first)/(samples[first+1].first-samples[first].first);
        float t2=t*t;
        float t3=t2*t;
        float h0=2.0f*t3-3.0f*t2+1.0f;
        float h1=-2.0f*t3+3.0f*t2;
        float h2=t3-2.0f*t2+t;
        float h3=t3-t2;
        return h0*samples[first].second+
            h1*samples[first+1].second+
            h2*tangents[first]+
            h3*tangents[first+1];
    }

    std::vector<std::pair<float,float> > samples;
    std::vector<float> tangents;
    bool even_sorted_sample;
    float step_size;
};

#endif // PBRT_CORE_SPLINES_H
