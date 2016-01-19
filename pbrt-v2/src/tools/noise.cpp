#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "../core/splines.h"
#include "../core/perlin.h"

using namespace std;

static void usage() {
    fprintf(stderr, "usage: spline <input> <output.txt>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    //const char *infile = NULL;

    if (argc < 2) usage();
    //infile = argv[1];
    outfile = argv[1];

    Perlin noise(0.25,1.0,3,10);
    ofstream output(outfile);
    for(size_t i=0;i<100;i++){
        float x=(float)i/100.0f;
        output<<x<<" "<<noise.Evaluate(x)<<endl;
    }
    output.close();

    return 0;
}
