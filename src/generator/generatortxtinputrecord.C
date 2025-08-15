/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "generatortxtinputrecord.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dictionary.h"
#include "range.h"
#include "scalarfunction.h"
#include "generatorerror.h"


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <ostream>
#include <sstream>

GeneratorTXTInputRecord::GeneratorTXTInputRecord() : tokenizer(), record()
{ }

GeneratorTXTInputRecord::GeneratorTXTInputRecord(const GeneratorTXTInputRecord &src) : tokenizer(),
    record(src.record), lineNumber(src.lineNumber)
{
    tokenizer.tokenizeLine(this->record);
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }
}

GeneratorTXTInputRecord::GeneratorTXTInputRecord(int linenumber, std::string source) : tokenizer(),
    record(std::move(source) ), lineNumber(linenumber)
{
    tokenizer.tokenizeLine(this->record);
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = false;
    }
}

GeneratorTXTInputRecord &
GeneratorTXTInputRecord::operator = ( const GeneratorTXTInputRecord &src )
{
    this->record = src.record;
    tokenizer.tokenizeLine(this->record);
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }

    return * this;
}

static inline void trim_right(std::string &s) {
    while ( !s.empty() && ( s.back() == '\r' || s.back() == '\n' ) ) {
        s.pop_back();
    }
}


void GeneratorTXTInputRecord::giveRawLine(std::string &dst) const
{
    dst = this->record;   // the original line as read by the reader
    trim_right(dst);      // drop trailing \r/\n so fopen gets a clean name
}

void
GeneratorTXTInputRecord::setRecordString(std::string newRec)
{
    this->record = std::move(newRec);
    tokenizer.tokenizeLine(this->record);
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = false;
    }
}

void
GeneratorTXTInputRecord::giveRecordKeywordField(std::string &answer, int &value)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std::string(tokenizer.giveToken(1) );
        setReadFlag(1);
        auto ptr = scanInteger(tokenizer.giveToken(2), value);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, "RecordID", lineNumber);
        }
        setReadFlag(2);
    } else {
        throw BadFormatInputException(*this, "RecordID", lineNumber);
    }
}


void
GeneratorTXTInputRecord::giveRecordKeywordField(std::string &keyword, std::string &value)
{
    if ( tokenizer.giveNumberOfTokens() >= 2 ) {
        keyword = std::string(tokenizer.giveToken(1) );
        value   = std::string(tokenizer.giveToken(2) );
        setReadFlag(1);
        setReadFlag(2);
    } else {
        throw BadFormatInputException(*this, "RecordID", lineNumber);
    }
}


void
GeneratorTXTInputRecord::giveRecordKeywordField(std::string &answer)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std::string(tokenizer.giveToken(1) );
        setReadFlag(1);
    } else {
        throw BadFormatInputException(*this, "RecordID", lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(int &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        auto ptr = scanInteger(tokenizer.giveToken(indx + 1), answer);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(double &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        auto ptr = scanDouble(tokenizer.giveToken(indx + 1), answer);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(bool &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int val;
        auto ptr = scanInteger(tokenizer.giveToken(indx + 1), val);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
        answer = val != 0;
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(std::string &answer, InputFieldType id)
{
    int indx = 0;
    if ( id ) {
        if ( ( indx = this->giveKeywordIndx(id) ) == 0 ) {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        indx++;
    } else {
        indx = 1;
    }

    const char *_token = tokenizer.giveToken(indx);
    if ( _token ) {
        answer = std::string(tokenizer.giveToken(indx) );
        setReadFlag(indx);
    } else {
        answer = "";
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(oofem::IntArray &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int size;
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            int value;
            ptr = scanInteger(tokenizer.giveToken(indx + i), value);
            if ( ptr == nullptr || * ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            answer.at(i) = value;
            setReadFlag(indx + i);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(oofem::FloatArray &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int size;
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            double value;
            auto ptr = scanDouble(tokenizer.giveToken(indx + i), value);
            if ( ptr == nullptr || * ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            answer.at(i) = value;
            setReadFlag(indx + i);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(oofem::FloatMatrix &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int nrows, ncols;
        setReadFlag(indx);

        auto ptr = scanInteger(tokenizer.giveToken(++indx), nrows);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        ptr = scanInteger(tokenizer.giveToken(++indx), ncols);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);

        if ( readMatrix(tokenizer.giveToken(++indx), nrows, ncols, answer) == 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(std::vector < std::string > & answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int size;
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }
        answer.reserve(size);
        setReadFlag(indx);
        for ( int i = 1; i <= size; i++ ) {
            answer.push_back(tokenizer.giveToken(indx + i) );
            setReadFlag(indx + i);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(oofem::Dictionary &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        int size;
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == nullptr || * ptr != 0 ) {
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);

        answer.clear();
        for ( int i = 1; i <= size; i++ ) {
            int key = 0;
            const char *token = tokenizer.giveToken(++indx);
            auto ptr1 = scanInteger(token, key);
            if ( ptr1 == nullptr || * ptr1 != 0 ) {
                key = token [ 0 ];
                // throw BadFormatInputException(*this, id, lineNumber);
            }
            double value;
            setReadFlag(indx);
            auto ptr = scanDouble(tokenizer.giveToken(++indx), value);
            if ( ptr == nullptr || * ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            answer.add(key, value);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(std::list < oofem::Range > & list, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        int li, hi;
        setReadFlag(indx);
        const char *rec = tokenizer.giveToken(++indx);
        if ( * rec != '{' ) {
            printf("missing left '{'");
            list.clear();
            throw BadFormatInputException(*this, id, lineNumber);
        }

        setReadFlag(indx);
        rec++;
        // read ranges
        while ( readRange(& rec, li, hi) ) {
            oofem::Range range(li, hi);
            list.push_back(range);
        }

        // skip whitespaces after last range
        while ( isspace(* rec) ) {
            rec++;
        }

        // test for enclosing bracket
        if ( * rec != '}' ) {
            printf("missing end '}'");
            list.clear();
            throw BadFormatInputException(*this, id, lineNumber);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

void
GeneratorTXTInputRecord::giveField(oofem::ScalarFunction &answer, InputFieldType id)
{
    const char *rec;
    int indx = this->giveKeywordIndx(id);

    if ( indx ) {
        setReadFlag(indx);
        rec = tokenizer.giveToken(++indx);
        if ( * rec == '@' ) {
            // reference to function
            int refVal;
            auto ptr = scanInteger(rec + 1, refVal);
            if ( ptr == nullptr || * ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }
            setReadFlag(indx);
            answer.setReference(refVal);
        } else if ( * rec == '$' ) {
            // simple expression
            std::string expr;

            expr = std::string(tokenizer.giveToken(indx) );
            setReadFlag(indx);
            std::string _v = expr.substr(1, expr.size() - 2);

            answer.setSimpleExpression(_v); // get rid of enclosing '$'
        } else {
            double val;
            auto ptr = scanDouble(tokenizer.giveToken(indx), val);
            if ( ptr == nullptr || * ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            answer.setValue(val);
        }
    } else {
        throw MissingKeywordInputException(*this, id, lineNumber);
    }
}

bool
GeneratorTXTInputRecord::hasField(InputFieldType id)
{
    //returns nonzero if id is present in source
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
    }

    return ( indx > 0 ) ? true : false;
}

void
GeneratorTXTInputRecord::printYourself()
{
    printf("%s", this->record.c_str() );
}

const char *
GeneratorTXTInputRecord::scanInteger(const char *source, int &value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == nullptr ) {
        value = 0;
        return nullptr;
    }

    value = strtol(source, & endptr, 10);
    return endptr;
}

const char *
GeneratorTXTInputRecord::scanDouble(const char *source, double &value)
{
    //
    // reads double value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == nullptr ) {
        value = 0;
        return nullptr;
    }

    value = strtod(source, & endptr);
    return endptr;
}

int
GeneratorTXTInputRecord::giveKeywordIndx(const char *kwd)
{
    int ntokens = tokenizer.giveNumberOfTokens();
    for ( int i = 1; i <= ntokens; i++ ) {
        if ( strcmp(kwd, tokenizer.giveToken(i) ) == 0 ) {
            return i;
        }
    }

    return 0;
}

void
GeneratorTXTInputRecord::finish(bool wrn)
{
    if ( !wrn ) {
        return;
    }

    std::ostringstream buff;
    bool pf = true, wf = false;
    int ntokens = tokenizer.giveNumberOfTokens();
    for ( int i = 0; i < ntokens; i++ ) {
        //fprintf (stderr, "[%s] ", tokenizer.giveToken(i+1));
        if ( !readFlag [ i ] ) {
            if ( pf ) {
                buff << "Unread token(s) detected in the following record\n\"";
                for ( int j = 0; j < 40; j++ ) {
                    if ( this->record [ j ] == '\n' || this->record [ j ] == '\0' ) {
                        break;
                    } else {
                        buff << this->record [ j ];
                    }
                }
                if ( this->record.size() > 41 ) {
                    buff << "...";
                }
                buff << "\":\n";

                pf = false;
                wf = true;
            }

            buff << "[" << tokenizer.giveToken(i + 1) << "]";
        }
    }

    if ( wf ) {
        printf(buff.str().c_str() );
    }
}

int
GeneratorTXTInputRecord::readRange(const char **helpSource, int &li, int &hi)
{
    char *endptr;
    // skip whitespaces
    while ( isspace(* * helpSource) ) {
        ( * helpSource )++;
    }

    // test if character is digit
    if ( isdigit(* * helpSource) ) {
        // digit character - read one value range
        li = hi = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        return 1;
    } else if ( * * helpSource == '(' ) {
        // range left parenthesis found
        ( * helpSource )++;
        // read lower index
        li = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        // test whitespaces
        if ( * * helpSource != ' ' && * * helpSource != '\t' ) {
            printf("unexpected token while reading range value");
            return 0;
        }

        // read end index
        hi = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        // skip whitespaces
        while ( isspace(* * helpSource) ) {
            ( * helpSource )++;
        }

        // test for enclosing bracket
        if ( * * helpSource == ')' ) {
            ( * helpSource )++;
            return 1;
        } else {
            printf("end ')' missing while parsing range value");
            return 0;
        }
    }

    return 0;
}

int
GeneratorTXTInputRecord::readMatrix(const char *helpSource, int r, int c, oofem::FloatMatrix &ans)
{
    const char *endptr = helpSource;

    if ( helpSource == NULL ) {
        ans.clear();
        return 0;
    }

    ans.resize(r, c);
    // skip whitespaces
    while ( isspace(* endptr) ) {
        ( endptr )++;
    }

    if ( * endptr == '{' ) {
        // range left parenthesis found
        ( endptr )++;
        // read row by row separated by semicolon
        for ( int i = 1; i <= r; i++ ) {
            for ( int j = 1; j <= c; j++ ) {
                endptr = scanDouble(endptr, ans.at(i, j) );
            }

            if ( i < r ) {
                // skip whitespaces
                while ( isspace(* endptr) ) {
                    ( endptr )++;
                }

                // test for row terminating semicolon
                if ( * endptr == ';' ) {
                    ( endptr )++;
                } else {
                    printf("missing row terminating semicolon");
                    return 0;
                }
            }
        }

        // skip whitespaces
        while ( isspace(* endptr) ) {
            ( endptr )++;
        }

        // test for enclosing bracket
        if ( * endptr == '}' ) {
            return 1;
        } else {
            printf("end '}' missing while parsing matrix value");
            return 0;
        }
    } else {
        return 0;
    }
}
