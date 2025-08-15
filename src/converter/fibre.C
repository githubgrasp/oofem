
#include "fibre.h"
#include "vertex.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "line.h"
#include <vector>


#endif


Fibre :: Fibre(int n, Grid *aGrid) : GridComponent(n, aGrid)//, coordinates()
{
    this->number = n;
    m_TOL=1.0e-10;
    

};

Fibre :: ~Fibre()
// Destructor.
{


};

IRResultType
Fibre :: initializeFrom(InputRecord *ir)
// Gets from the source Fibre from the data file all the data of the receiver.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    
    IR_GIVE_FIELD(ir, m_endpoints, IFT_Fibre_endpoints, "endpoints");
    IR_GIVE_FIELD(ir, diameter, IFT_Fibre_diameter, "diameter");
    oofem::FloatArray tool(3);
    coordP=tool;
    coordQ=tool;// to allocate the size
    direction_vector=tool;
    
    coordP.at(1)=m_endpoints.at(1);
    coordP.at(2)=m_endpoints.at(2);
    coordP.at(3)=m_endpoints.at(3);
    coordQ.at(1)=m_endpoints.at(4);
    coordQ.at(2)=m_endpoints.at(5);
    coordQ.at(3)=m_endpoints.at(6);



    PointP =  grid->createInterNode(coordP);
    PointQ =  grid->createInterNode(coordQ);
    
    //computation of the direction vector
    length=computedistance(coordP,coordQ);
    direction_vector=(1./length)*(coordQ-coordP);
    return IRRT_OK;
    
};

// tool functions
void Fibre :: findDelaunayVertex()
{
    int N (m_IntersectionPoints.size());
    oofem::FloatArray coordsOne, coordsTwo;
    IntArray nodeCandidates;
    int index_SearchedPoint=1, index_SearchedPointTwo=1;
    if (N>0)
    {
        Vertex *interP;

        interP=grid->giveInterNode(m_IntersectionPoints.at(N-1));
        interP->giveCoordinates(coordsOne);
        // look for the nearest two Delaunay vertex
      
        // -> preselection of candidates nodes
        nodeCandidates=grid->findDelaunayNodesWithinBox(coordsOne, 2.*grid->giveDiameter());
        
        int NbDelaunayVertices(nodeCandidates.giveSize());
        if(NbDelaunayVertices<1 && grid->giveMeshType() == 0){//Debug
	  printf("error in descritization of fibre, preselection in localization of Delaunay first vertex must be wider");
	  exit(1);
        }
	else if(NbDelaunayVertices<1 && grid->giveMeshType() == 1){//Debug
	  printf("error in descritization of fibre, preselection in localization of Delaunay first vertex must be wider");
	  exit(1);
        }

        //initialisation
        double distance1 = 2*grid->giveDiameter();
	double distance;
        
        for (int i=1;i<=NbDelaunayVertices;i++){  
	  grid->giveDelaunayVertex(nodeCandidates.at(i))->giveCoordinates(coordsTwo);
	  distance= computedistance(coordsOne,coordsTwo);
            
	  if(distance<distance1){
	    index_SearchedPoint=nodeCandidates.at(i);
	    distance1=distance;
	  }
        }
	
        m_DelaunayVertices.push_back(index_SearchedPoint);
        
        }
    else {
      printf("error in findDelaunayVertex in discretization of fibres");
    }
    return;
}



void Fibre :: findIntersect()
{
  int NbIntPoints = m_IntersectionPoints.size();
  
  if (NbIntPoints==0){//First point

    //Use start point of fibre as first intersection
    m_IntersectionPoints.push_back(this->PointP->giveNumber());  

  }
  else {
      
        int NbDelVertices = m_DelaunayVertices.size();
        int num = m_DelaunayVertices.at(NbDelVertices-1);
        
        // finding of neighbourgh nodes
        IntArray NeighbLines;// warning : index start 0
        grid->giveDelaunayVertex(m_DelaunayVertices.at(NbDelVertices-1))->giveLocalLines(NeighbLines);
        
        IntArray nodes(2);
        int NbofNeighbour = NeighbLines.giveSize();
        IntArray Nb_neighbour(NbofNeighbour);

        for(int i=0;i<NbofNeighbour;i++)
        {
            grid->giveDelaunayLine(NeighbLines(i))->giveLocalVertices(nodes);

	    // we keep the number different from those of the initial node
            if (nodes.at(1)==num)
            {
                Nb_neighbour.at(i+1)=nodes.at(2);
            }
            else if (nodes.at(2)==num)
            {
                Nb_neighbour.at(i+1)=nodes.at(1);    
            }
            else
            {
                printf("error in the discretization of fibre, bad use of data");
            }
        }
        
        
        oofem::FloatArray coordJ,coordI,coordM,coordU;
        oofem::FloatArray diffIJ,diffMU,diffQU, coordS;
        double scalar1,scalar2;
        double rho_min(2),rho_candidate;//init at 2 car rho_min will necessarily be found smaller than 1
        int index_neighbour;
        
        grid->giveDelaunayVertex(m_DelaunayVertices.at(NbDelVertices-1))->giveCoordinates(coordI);
    
        grid->giveInterNode(m_IntersectionPoints.at(NbIntPoints-1))->giveCoordinates(coordU);
        

        //computation of rho for each of them
        for(int i=0;i<NbofNeighbour;i++)
        {
            
            grid->giveDelaunayVertex(Nb_neighbour.at(i+1))->giveCoordinates(coordJ);
            coordM=0.5*(coordI+coordJ);
            diffMU=coordM-coordU;
            diffIJ=coordI-coordJ;
            diffQU=coordQ-coordU;
            scalar1=dotProduct(diffIJ,diffMU,3);
            scalar2=dotProduct(diffIJ,diffQU,3);
            
            
            if (!scalar2==0) // either parallel lines...
            {
                rho_candidate=scalar1/scalar2;
                if (rho_candidate<rho_min && rho_candidate>0.+m_TOL)
                {
                    rho_min=rho_candidate;
                    index_neighbour=Nb_neighbour.at(i+1);
                }
            }
        }
        
        // we choose the smaller rho : thus we deduce the following cell to consider (neighbour)
        // and we create the intersection point with the facet
        

        if (rho_min<1)
            
        {   coordS=coordU+rho_min*diffQU;
            
            m_IntersectionPoints.push_back(grid->createInterNode(coordS)->giveNumber());
            
            m_DelaunayVertices.push_back(grid->giveDelaunayVertex(index_neighbour)->giveNumber());

        }
        else if (rho_min>=1) // here we have reached the end of the fibre, so we switch for the endpoint
        {
            m_IntersectionPoints.push_back(PointQ->giveNumber());
        }
        else
        {
            printf("error in the discretization of fibre");
        }
        
    }
    
    
    return;
    
}

void Fibre :: placeReinforcementNode()
{
    int NbReiPoints = m_ReinforcementPoints.size();
    
    
    // coordinates of the reinforcement node
    oofem::FloatArray coordR,coord1,coord2;
    
    int NbIntPoints = m_IntersectionPoints.size();
    grid->giveInterNode(m_IntersectionPoints.at(NbIntPoints-2))->giveCoordinates(coord1);
    grid->giveInterNode(m_IntersectionPoints.at(NbIntPoints-1))->giveCoordinates(coord2);
    
    coordR=0.5*(coord1+coord2);

    
    m_ReinforcementPoints.push_back(grid->createReinfNode(coordR)->giveNumber() );
    
    // computing of the associated L_end (length to the nearest endpoint
    double cL;
    cL=computedistance(coordR,coordP);
    m_listOfL_end.push_back(std::min(cL,length-cL));
                            
    //control
    //coordR.printYourself();
    //printf("point of reinf placed \n ----------- \n ");
    //int entier;
    //std::cin>>entier;

    return;
}

// main discretisation function
void Fibre :: discretizeYouself()
{
    int task(1);
    this->findIntersect(); // first call to place first endpoint in intersectionPoints list
    this->findDelaunayVertex();
    
    int i(0);
    
    while(task==1){
      i++;
      if (i>1000){
	printf("Error in discretization : max number of iteration reached.");
	printf("Consider to adjust parameter m_TOL to avoid instabilities");
	exit(1);
      }
      
      int NbIntPoints(m_IntersectionPoints.size());
      if (m_IntersectionPoints.at(NbIntPoints-1)==PointQ->giveNumber()){
	task=0;
      }
      else{
	this->findIntersect();
	this->placeReinforcementNode();
      }
    }
    
    return;
}


Fibre *Fibre :: ofType()
// Returns a new DofManager, which has the same number than the receiver,
// but belongs to aClass (Node, ElementSide,..).
{
    Fibre *fibre;
    
    fibre = new Fibre(number,grid);
    
    return fibre;
}

//accessors
int Fibre ::giveNumberReinforcementNode(int i)
{
    return m_ReinforcementPoints.at(i-1);
}

int Fibre ::NbOfReinfNodes()
{
    return m_ReinforcementPoints.size();
}

int Fibre ::giveNumberDelaunayNode(int i)
{
    return m_DelaunayVertices.at(i-1);

}

int Fibre ::NbOfDelNodes()
{
    return m_DelaunayVertices.size();

}

int Fibre ::giveNumberIntersectionPoint(int i)
{
    return m_IntersectionPoints.at(i-1);
    
}

int Fibre ::NbOfIntersectionPoints()
{
    return m_IntersectionPoints.size();
    
}

double Fibre ::giveL_end(int i)
{
    return m_listOfL_end.at(i-1);
}

int Fibre ::NbOfL_end()
{
    return m_listOfL_end.size();
}

// static functions

Vertex *Fibre::reinforcementLocalizer(oofem::FloatArray coord, Grid* agrid, double TOL)
{  // find the reinforcement node located in vicinity of the specified coordinates (with TOL as tolerance)
    // error if no node found
    // TOL must be chosen enough small, or the returned node will be the first found
    bool test(false);
    Vertex *SearchedNode;
    oofem::FloatArray coordN;
    double distance;
    
    
    for(int i=0;i<agrid->giveNumberOfReinforcementNode();i++)
    {
        agrid->giveReinforcementNode(i+1)->giveCoordinates(coordN);
        
        distance=computedistance(coord,coordN);
        if (distance<=TOL)
        {
            SearchedNode=agrid->giveReinforcementNode(i+1);
            test=true;
            break;
        }
        
    }
    
    
    if (test==false){
      
      // search for the min distance        
      double distance2;
      distance=100;
        
      for(int i=0;i<agrid->giveNumberOfReinforcementNode();i++)
        {
	  agrid->giveReinforcementNode(i+1)->giveCoordinates(coordN);
            
            distance2=computedistance(coord,coordN);
            if (distance2<=distance)
            {
	      distance=distance2;
	      SearchedNode=agrid->giveReinforcementNode(i+1);	      
            }
            
        }
        
        
        printf("\n reinforcementLocalizer : no node found. Is the input really periodic? for info : Tolerance is %e \n", TOL);
        printf("\n for info, the nearest node was found at %e \n", distance);
	//    exit(1);
    }
    
    return SearchedNode;
}


double Fibre::computedistance(oofem::FloatArray coordsOne,oofem::FloatArray coordsTwo)
{
    double distance;
    distance=sqrt( pow(coordsOne.at(1) - coordsTwo.at(1), 2.) +
                   pow(coordsOne.at(2) - coordsTwo.at(2), 2.) +
                   pow(coordsOne.at(3) - coordsTwo.at(3), 2.) );
    return distance;
}

