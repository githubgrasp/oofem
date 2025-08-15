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

#include <string>
#include <algorithm>
#include <cctype>
#include "generatorinputrecord.h"

void
ConverterInputRecord::giveOptionalField(int &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(double &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(bool &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(std::string &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(oofem::FloatArray &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(oofem::IntArray &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(oofem::FloatMatrix &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(std::vector < std::string > & answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(oofem::Dictionary &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(std::list < oofem::Range > & answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}

void
ConverterInputRecord::giveOptionalField(oofem::ScalarFunction &answer, InputFieldType id)
{
    if ( this->hasField(id) ) {
        try {
            this->giveField(answer, id);
        } catch(MissingKeywordInputException &) { }
    }
}


InputException::InputException(const ConverterInputRecord &ir, std::string keyword, int number) :
    record(ir.giveRecordAsString() ), keyword(std::move(keyword) ), number(number)
{ }


MissingKeywordInputException::MissingKeywordInputException(const ConverterInputRecord &ir, std::string kw, int n) :
    InputException(ir, std::move(kw), n)
{
    msg = "Missing keyword \"" + keyword + "\" on input " + std::to_string(number) + \
          "\nRecord: \"" + record.substr(0, 50) + ( record.size() > 50?"...":"" ) + "\"";
}


BadFormatInputException::BadFormatInputException(const ConverterInputRecord &ir, std::string kw, int n) :
    InputException(ir, std::move(kw), n)
{
    msg = "Bad format for keyword \"" + keyword + "\" on input " + std::to_string(number) + \
          "\nRecord: \"" + record.substr(0, 50) + ( record.size() > 50?"...":"" ) + "\"";
}


ValueInputException::ValueInputException(const ConverterInputRecord &ir, std::string kw, const std::string &reason) :
    InputException(ir, std::move(kw), -1)
{
    msg = "Value input error for keyword \"" + keyword + "\"" + \
          "\nReason: \"" + reason + "\"\n" + \
          "\nRecord: \"" + record.substr(0, 50) + ( record.size() > 50?"...":"" ) + "\"";
}

void ConverterInputRecord::giveRawLine(std::string &out) const
{
    out = this->rawLine;  // assuming you have stored the raw line somewhere

    // --- trim leading whitespace
    out.erase(out.begin(),
              std::find_if(out.begin(), out.end(),
                           [] ( unsigned char ch ) { return !std::isspace(ch);
                           }) );

    // --- trim trailing whitespace
    out.erase(std::find_if(out.rbegin(), out.rend(),
                           [] ( unsigned char ch ) { return !std::isspace(ch);
                           }).base(),
              out.end() );
}


const char * MissingKeywordInputException::what() const noexcept
{
    return msg.c_str();
}


const char * BadFormatInputException::what() const noexcept
{
    return msg.c_str();
}

const char * ValueInputException::what() const noexcept
{
    return msg.c_str();
}
