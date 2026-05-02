#include "rect.h"


Rect::Rect(int n, Grid *aGrid) : Region(n, aGrid)
{
    this->number = n;
}


void Rect::defineBoundaries(oofem::FloatArray &boundaries)
{
    boundaries.resize(6);
    boundaries.at(1) = box.at(1);   // xmin
    boundaries.at(2) = box.at(3);   // xmax
    boundaries.at(3) = box.at(2);   // ymin
    boundaries.at(4) = box.at(4);   // ymax
    boundaries.at(5) = 0.;          // zmin
    boundaries.at(6) = 0.;          // zmax
}


bool Rect::contains(double x, double y, double tol) const
{
    return x > box.at(1) - tol && x < box.at(3) + tol &&
           y > box.at(2) - tol && y < box.at(4) + tol;
}
