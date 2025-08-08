#ifndef refineprism_h
#define refineprism_h

#include "refinement.h"

#include "grid.h"

#include "floatarray.h"
#include "intarray.h"

#include "inputrecord.h"

#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


#define _IFT_RefinePrism_box "box"
#define _IFT_RefinePrism_refine "refine"


class RefinePrism : public Refinement
{

protected:
    /// Array storing nodal coordinates.
  int number;
  double refinement;
  oofem::FloatArray box;

  
public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    RefinePrism(int n, Grid* aGrid);                      // constructor
    /// Destructor.
    ~RefinePrism();                                           // destructor

  virtual double giveDiameter(oofem::FloatArray &coord) override;
    
    RefinePrism* ofType();

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "RefinePrism"; }

  void initializeFrom(GeneratorInputRecord &ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();

};

#endif // refineprism_h






