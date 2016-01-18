#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <core/splines.h>

using namespace std;

static void usage() {
    fprintf(stderr, "usage: spline <input> <output.txt>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    const char *infile = NULL;

    if (argc < 3) usage();
    infile = argv[1];
    outfile = argv[2];

    Catmull_Rom spline;
    ifstream input(infile);
    float x,y;
    while(input>>x){
        input>>y;
        spline.Add_Sample(x,y);
    }
    input.close();
    spline.Build();
    ofstream output(outfile);
    float start_x=spline.samples[0].first;
    size_t output_sample_number=16;
    float step=(spline.samples[spline.samples.size()-1].first-spline.samples[0].first)/(float)(output_sample_number-1);
    for(size_t i=0;i<output_sample_number;i++){
        float x=start_x+step*(float)i;
        output<<x<<" "<<spline.Evaluate(x)<<endl;
    }
    output.close();

    return 0;
}
