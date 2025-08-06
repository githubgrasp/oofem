#ifndef octreegridlocalizer_h
#define octreegridlocalizer_h

#include "gridlocalizer.h"
//#include "compiler.h"

#ifndef __MAKEDEPEND
 #include <set>
 #include <list>
#endif

class Grid;
class Element;
class TimeStep;
class OctreeGridLocalizer;
/// Max desired number of nodes per octant
#define OCTREE_MAX_NODES_LIMIT 50
/// Max octree depth
#define OCTREE_MAX_DEPTH 15


/**
 * Class representing the octant of octree.
 * It maintains the link to parent cell or if it it the root cell, this link pointer is set to NULL.
 * Maintains links to possible child octree cells as well as its position and size.
 * Also list of node numbers contained in given octree cell can be maintained if cell is terminal cell.
 *
 */
class oofemOctantRec
{
protected:

    /// link to octree class
    OctreeGridLocalizer *localizer;
    /// link to parent cell record
    oofemOctantRec *parent;
    /// link to octant childs
    oofemOctantRec *child [ 2 ] [ 2 ] [ 2 ];
    /// octant origin coordinates
    FloatArray origin;
    /// octant size
    double size;

    /// octant node list
    std :: list< int > *nodeList;
    /// element list, containing all elements having ip in cell
    std :: set< int > *elementList;

public:
    enum boundingBoxStatus { BBS_OUTSIDECELL, BBS_INSIDECELL, BBS_CONTAINSCELL };

    /// constructor
    oofemOctantRec(OctreeGridLocalizer *loc, oofemOctantRec *parent, FloatArray &origin, double size);
    /// destructor
    ~oofemOctantRec();

    /// returns reference to parent; NULL if root
    oofemOctantRec *giveParent() { return this->parent; }
    /// returns the cell origin
    void giveOrigin(FloatArray &answer) { answer = this->origin; }
    /// returns cell size
    double giveSize() { return this->size; }
    /** Returns nonzero if octant contains given point.
     * If not 3-coordinates are given, then missing coordinates are
     * not included in test */
    int containsPoint(const FloatArray &coords);
    /// returns the child cell with given local cell coordinates
    oofemOctantRec *giveChild(int xi, int yi, int zi);
    /**
     * Returns the child containing given point.
     * If not full 3d coordinates are provided, then only provided coordinates are taken into account,
     * assuming remaining to be same as origin.
     * If point is contained by receiver, coressponding child is set, otherwise
     * child is set to NULL.
     * @return 1 if o.k, -1 if no childs exists (receiver is terminal octant), -2 point out of receiver volume
     */
    int giveChildContainingPoint(oofemOctantRec **child, const FloatArray &coords);
    /// Returns nonzero if octant is terminal one (no children)
    int isTerminalOctant();
    /// Return reference to node List
    std :: list< int > *giveNodeList();
    /// Return reference to IPelement set
    std :: set< int > *giveIPElementList();

    /**
     * Divide receiver further, creating corresponding childs
     */
    int divideLocally(int level, const IntArray &octantMask);
    /**
     * Test if receiver within bounding box (sphere)
     * @returns boundingBoxStatus status
     */
    boundingBoxStatus testBoundingBox(const FloatArray &coords, double radius);
    /**
     * Adds given element to cell list of elements having IP within this cell.
     */
    void addElementIP(int elementNum) { this->giveIPElementList()->insert(elementNum); }
    /**
     * Adds given Node to node list of nodes contained by receiver
     */
    void addNode(int nodeNum) { this->giveNodeList()->push_back(nodeNum); }
    /**
     * Clears and deletes the nodeList
     */
    void deleteNodeList() {
        if ( elementList ) { delete elementList; }

        elementList = NULL;
    }
};





/**
 * The implementation of grid localizer based on octree technique.
 * The basic task is to provide grid information and localization for grid, to which receiver is associated.
 * The octree data structure is build over grid nodes, the element type seraches are completed using
 * nodal conectivity informations provided by ConTable.
 * Typical services include searching the closes node to give position, serching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class OctreeGridLocalizer : public GridLocalizer
{
protected:
    /// Root cell of octree
    oofemOctantRec *rootCell;
    /// Octree degenerate mask
    IntArray octreeMask;
    /// Flag inficating elementIP tables are initialized
    int elementIPListsInitialized;
public:
    /// Constructor
    OctreeGridLocalizer(int n, Grid *d) : GridLocalizer(n, d), octreeMask(3) {
        rootCell = NULL;
        elementIPListsInitialized = 0;
    }
    /// Destructor - deletes the octree tree
    ~OctreeGridLocalizer() { if ( rootCell ) { delete rootCell; } }


    /**
     * Initialize receiver data structure if not done previously
     * Current implementation calls and returns the buildOctreeDataStructure service response.
     */
    int init(bool force = false);


    int checkNodesWithinBox(const FloatArray &coords, const double radius);

    /**
     * Returns container (list) of all grid nodes within given box.
     * @param NODESet answer containing the list of nodes meeting the criteria
     * @param coords center of box of interest
     * @param radius radius of bounding sphere
     */
    void giveAllNodesWithinBox(nodeContainerType &nodeList, const FloatArray &coords, const double radius);



    /**
     * Returns container (list) of nodes within given box and given root cell.
     * @param elemSet answer containing the list of nodes meeting the criteria
     * @param currentCell the starting cell to be transversed
     * @param coords center of box of interest
     * @param radius radius of bounding sphere
     */
    void giveNodesWithinBox(nodeContainerType &nodeList, oofemOctantRec *currentCell,
                            const FloatArray &coords, const double radius);


    void insertSequentialNode(int nodeNum, const FloatArray &coords);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "OctreeGridLocalizer"; }
    /** Returns the octreeMask value given by the index */
    int giveOctreeMaskValue(int indx) { return octreeMask.at(indx); }


protected:
    /**
     * Buids the underlying octree data structure.
     * The desired tree level is determined by following rules:
     * - if number of nodes exceed treshold, cell is subdivided
     * - there is maximal octree level.
     * - In current implementation, the neighbour cell size diference is allowed to be > 2
     */
    int buildOctreeDataStructure();

    /**
     * Inserts the given node (identified by its number and position) to the octree structure.
     * The tree is traversed until terminal octant containing given position is found and node is then inserted
     * into octant nodal list. If there is too mach nodes per cell, this is subdivided further and
     * its assigned nodes are propagated to childs.
     * @param rootCell starting cell for insertion
     * @param nodeNum node number
     * @param coords corresponding node coordinates
     * @return nonzero if node insertion was succesfull
     */
    int insertNodeIntoOctree(oofemOctantRec *rootCell, int nodeNum, const FloatArray &coords);

    /**
     * Finds the terminal octant containing the given point.
     * @param startCell cell used to start search
     * @param coords coordinates of point of interest
     * @return pointer to terminal octant, NULL if point outside startingCell
     */
    oofemOctantRec *findTerminalContaining(oofemOctantRec *startCell, const FloatArray &coords);


    /**
     * Returns the element close to given point as element which center is closest to the given point.
     * Only elements with IP in given cell are considered.
     * The search is done only for given cell and its childs, skipping the given child from search
     * @param cell top level cell to search
     * @param coords point coordinates
     * @param minDist distance from the center of returned element
     * @param regionList only elements within given regions are considered, if NULL all regions are considered.
     */
    void giveElementCloseToPointWithinOctant(oofemOctantRec *cell, const FloatArray &coords,
                                             double &minDist, Element **answer,
                                             const IntArray *regionList);
    /**
     * Returns the octree depth for given cell. The depth is not parameter of octree cells, but is
     * computed from cell size and root cell size.
     */
    int giveCellDepth(oofemOctantRec *cell);
    /**
     * Determines the max tree depth computed for given tree cell and its childs.
     * To obtain ttotal tree depth, root cell should be supplied.
     * The tree depth is always measured from the root cell.
     * Before call, maxDepth should be set to zero.
     */
    void giveMaxTreeDepthFrom(oofemOctantRec *root, int &maxDepth);
    /**
     * Builds the list of terminal cells contained within given box (coords, radius), starting from given currentCell.
     * @param cellList list of terminal cell pointers contained by bbox
     * @param coords center of box of interest
     * @param radius radius of bounding sphere
     * @param currentCell starting cell
     */
    void giveListOfTerminalCellsInBoundingBox(std :: list< oofemOctantRec * > &cellList, const FloatArray &coords,
                                              const double radius, oofemOctantRec *currentCell);
};

#endif // octreelocalizer_h
