#include "refineprism.h"
#include "curve.h"
#include "vertex.h"
#include "generatorerror.h"


#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdlib.h>
 #include "octreegridlocalizer.h"
#endif


RefinePrism::RefinePrism(int n, Grid *aGrid) : Refinement(n, aGrid) //, coordinates()
{
}

RefinePrism::~RefinePrism()
{}

void RefinePrism::initializeFromTokens(std::istringstream &iss)
{
    refinement = 1.;

    bool gotBox = false;
    std::string tok;
    while ( iss >> tok ) {
        if ( tok == "box" ) {
            int n;
            iss >> n;
            box.resize(n);
            for ( int i = 1; i <= n; ++i ) {
                iss >> box.at(i);
            }
            gotBox = true;
        } else if ( tok == "refine" ) {
            iss >> refinement;
        } else {
            generator::errorf("RefinePrism::initializeFromTokens: unknown keyword '%s'", tok.c_str());
        }
    }
    if ( !gotBox ) {
        generator::error("RefinePrism::initializeFromTokens: missing 'box' keyword");
    }
}

double
RefinePrism::giveDiameter(oofem::FloatArray &coord) {
    //check of coord is in prism
    double diam = 0.;

    if ( coord.at(1) >= box.at(1) && coord.at(2) >= box.at(2) && coord.at(3) >= box.at(3) && coord.at(1) <= box.at(4) && coord.at(2) <= box.at(5) && coord.at(3) <= box.at(6) ) {//inside the box
        diam = this->refinement * grid->diameter;
    } else   {
        diam = grid->diameter;
    }

    return diam;
}
