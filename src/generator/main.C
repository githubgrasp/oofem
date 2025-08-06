#include "logger.h"

/* Default oofem loggers */
Logger oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include "flotarry.h"
#include "intarray.h"
#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"
#include "vertex.h"
#include "grid.h"


/// Returns class name of the receiver.
 char *giveClassName() { return "Main"; }

int main(int argc, char* argv[]){
  if(argc != 2){
    printf("wrong number of arguments: name of input file is needed\n");
    std::exit(1);
  }
  
  Grid *grid = new Grid(1);

  //These two things are required by the macro
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                // Required by IR_GIVE_FIELD macro
  
  char inputFileName [ MAX_FILENAME_LENGTH + 10 ], buff [ MAX_FILENAME_LENGTH ];
  char dataOutputFileName [ MAX_FILENAME_LENGTH + 10 ];
  strcpy(inputFileName, argv[1]);
  OOFEMTXTDataReader dr(inputFileName);

  //Read the output file name
  InputRecord *ir = dr.giveInputRecord(DataReader :: IR_outFileRec, 1);
  //  __keyword = NULL;
  result = ir->giveField(dataOutputFileName, MAX_FILENAME_LENGTH, IFT_outfile, NULL);
  if ( result != IRRT_OK ) {
    IR_IOERR("", __proc, IFT_outfile, "Output file record", ir, result);
  }
 
  FILE *outputStream;  
  if ( ( outputStream = fopen(dataOutputFileName, "w") ) == NULL ) {
    printf("Can't open output file %s", dataOutputFileName);
    exit(1);
  }
  
  grid->instanciateYourself(&dr);
  
  grid->generatePoints();
  
  grid->giveOutput(outputStream);

  printf("numberOfVertices = %d\n", grid->giveNumberOfVertices());
  printf("Point generation completed\n");
  return 1;
}

