
#ifndef inputrecord_h
#define inputrecord_h


#include "intarray.h"
#include "floatarray.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
#endif

/**
 * Type defining the return values of InputRecord reading operations.
 * IRRT_OK the corresponding value to given keyword was successfully read.
 *       the answer parameter contains the value.
 * IRRT_NOTFOUND the keyword is not found; the answer is not modified
 * IRRT_BAD_FORMAT the keyword was found but the record is not correctly formated.
 */
enum IRResultType { IRRT_OK = 0, IRRT_NOTFOUND, IRRT_BAD_FORMAT };

/**
 * Enumeration type used to determine particular field in record
 */
enum InputFieldType {
    IFT_RecordIDField,
    IFT_test,
    IFT_EngngModel_probdescription,
    IFT_EngngModel_problem,
    IFT_Vertex_coords,
    IFT_Vertex_refine,
    IFT_Vertex_radius,
    IFT_Vertex_randomswitch,
    IFT_outfile,
    IFT_type,
    IFT_diam,
    IFT_dens,
    IFT_maxiter,
    IFT_nvertex,
    IFT_ncurve,
    IFT_nsurface,
    IFT_nregion,
    IFT_regadditor,
    IFT_ranint,
    IFT_Curve_vertices,
    IFT_Curve_refine,
    IFT_Surface_curves,
    IFT_Surface_refine,
    IFT_Region_surfaces,
    IFT_Region_refine,
    IFT_Sphere_diameter,
    IFT_Sphere_refine,
    IFT_RefinePrism_box
};


/**
 * Macro simplifying the erorr reporting.
 */
#define IR_IOERR(__class, __proc, __id, __keyword, __ir, __result)      \
        __ir->report_error(__class, __proc, __id, __keyword, __result, __FILE__, __LINE__);

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the erorr reporting.
 */
#define IR_GIVE_FIELD(__ir, __value, __id, __kwd) result = __ir->giveField(__value, __id, __kwd); \
        if ( result != IRRT_OK ) { IR_IOERR(giveClassName(), __proc, __id, __kwd, __ir, result); }

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the optional
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the erorr reporting.
 */
#define IR_GIVE_OPTIONAL_FIELD(__ir, __value, __id, __kwd) result = __ir->giveOptionalField(__value, __id, __kwd); \
        if ( result != IRRT_OK ) { IR_IOERR(giveClassName(), __proc, __id, __kwd, __ir, result); }

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the erorr reporting.
 */
#define IR_GIVE_FIELD2(__ir, __value, __id, __kwd, __opt) result = __ir->giveField(__value, __opt, __id, __kwd); \
        if ( result != IRRT_OK ) { IR_IOERR(giveClassName(), __proc, __id, __kwd, __ir, result); }

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the optional
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the erorr reporting.
 */
#define IR_GIVE_OPTIONAL_FIELD2(__ir, __value, __id, __kwd, __opt)      \
        result = __ir->giveOptionalField(__value, __opt, __id, __kwd);        \
        if ( result != IRRT_OK ) { IR_IOERR(giveClassName(), __proc, __id, __kwd, __ir, result); }
/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory record keyword (__kwd)
 * and its number (__value param). Includes also the erorr reporting.
 */
#define IR_GIVE_RECORD_KEYWORD_FIELD(__ir, __name, __value, __opt)      \
        result = __ir->giveRecordKeywordField(__name, __value, __opt);        \
        if ( result != IRRT_OK ) { IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "RecordIDField", __ir, result); }



/**
 * Class representing the general Input Record. The input record consis of several fields.
 * Provides several requesting functions for reading field values. The derived classes of
 * Input record can represent database records or text file records, allowing the transparent
 * input operations.
 * The input record after init phase should "contain" all relevant data, so the input record should
 * resolve all dependencies. THis allows to create a copy of input record instance for later use
 * without the need to re-open input files (used for metasteps).
 */
class InputRecord
{
protected:
public:
    /// Constructor. Creates an empty input record.
    InputRecord();
    /// Copy constructor
    InputRecord(const InputRecord &);
    /// Destructor
    virtual ~InputRecord() { }
    // Assingnment operator
    InputRecord &operator=(const InputRecord &);

    /** Creates a newly allocated copy of the receiver */
    virtual InputRecord *GiveCopy() = 0;

    /**@name Compulsory field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the record id field  (type of record) and its corresponding number
    virtual IRResultType giveRecordKeywordField(char *answer, int &value, int maxchar) = 0;
    /// Reads the record id field  (type of record)
    virtual IRResultType giveRecordKeywordField(char *answer, int maxchar) = 0;
    /// Reads the integer field value
    virtual IRResultType giveField(int &answer, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the double field value
    virtual IRResultType giveField(double &answer, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the char* field value
    virtual IRResultType giveField(char *answer, int maxchar, const InputFieldType fieldI, const char *idString) = 0;
    /// Reads the FloatArray field value
    virtual IRResultType giveField(FloatArray &answer, const InputFieldType fieldI, const char *idString) = 0;
    /// Reads the IntArray field value
    virtual IRResultType giveField(IntArray &answer, const InputFieldType fieldID, const char *idString) = 0;

    /**@name Optional field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the integer field value
    virtual IRResultType giveOptionalField(int &answer, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the double field value
    virtual IRResultType giveOptionalField(double &answer, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the char* field value
    virtual IRResultType giveOptionalField(char *answer, int maxchar, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the FloatArray field value
    virtual IRResultType giveOptionalField(FloatArray &answer, const InputFieldType fieldID, const char *idString) = 0;
    /// Reads the IntArray field value
    virtual IRResultType giveOptionalField(IntArray &answer, const InputFieldType fieldID, const char *idString) = 0;

    /// Returns true if record contains field identified by idString keyword
    virtual bool         hasField(const InputFieldType fieldID, const char *idString) = 0;

    /// Returns error string corresponding to given value of IRResultType type
    const char *strerror(IRResultType);

    /// Prints the error message
    void report_error(const char *_class, const char *proc, const InputFieldType fieldID, const char *kwd,
                      IRResultType result, const char *file, int line);

    /** terminates the current record session and if flag is true warnin is printed for unscanned tokens */
    virtual void finish(bool wrn = true) = 0;
};

#endif // inputrecord_h
