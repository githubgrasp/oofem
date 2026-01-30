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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#ifndef mpmevaluator_h
#define mpmevaluator_h


#include <iostream>
#include <vector>
#include <string>
#include <stack>
#include <map>
#include <memory>
#include <functional>
#include <regex>
#include <sstream>
#include "floatmatrix.h"

using namespace std;

namespace oofem {

/**
 * --- PART 1: TYPE SYSTEM & METADATA ---
 */

enum class DataType { SCALAR, MATRIX, FIELD_INFO };

// Holds information about physical fields (e.g., displacement u)
struct FieldMetadata {
    std::string name;
    int dofs; 
};

// Represents the shape of an expression node during compilation
struct TypedShape {
    DataType type;
    int rows, cols;
    std::shared_ptr<FieldMetadata> field_meta;
    
    // Static builders for clear intent
    static TypedShape Scalar() { return {DataType::SCALAR, 1, 1, nullptr}; }
    static TypedShape Matrix(int r, int c) { return {DataType::MATRIX, r, c, nullptr}; }
    static TypedShape Field(std::shared_ptr<FieldMetadata> m) { 
        return {DataType::FIELD_INFO, m->dofs, 1, m}; 
    }
};

/**
 * --- PART 2: VIRTUAL MACHINE INSTRUCTIONS ---
 */

enum class OpCode { CALL_FUNC, MAT_MUL };

struct Instruction {
    OpCode op;
    int func_id;             // Index for external C++ providers
    std::vector<int> inputs; // Indices into the memory pool
    
    int lhs_idx = -1, rhs_idx = -1;        
    bool lhs_trans = false, rhs_trans = false; // Lazy transpose flags
    int output_idx;          
};

/**
 * --- PART 3: THE COMPILER ---
 * Handles logic for RPN conversion, literal parsing, and dimension validation.
 */

class MPMCompiler {
    struct FunctionSig {
        DataType returnType;
        std::vector<DataType> argTypes;
        std::function<TypedShape(const std::vector<TypedShape>&)> resolver;
        int provider_id;
    };

    std::map<std::string, int> precedence = {{"*", 2}, {"+", 1}, {".T", 4}, {".t", 4}};
    std::map<std::string, FunctionSig> registry;
    std::map<std::string, TypedShape> constants;

    // Checks if a string token is a valid number
    bool is_number(const std::string& s) {
        if (s.empty()) return false;
        char* endPtr = nullptr;
        std::strtod(s.c_str(), &endPtr);
        return endPtr != s.c_str() && *endPtr == '\0';
    }

    // Parses literals like "[1,2;3,4]" into Eigen matrices
    FloatMatrix parse_matrix_literal(std::string s) {
        s = std::regex_replace(s, std::regex(R"([\[\]])"), "");
        std::vector<std::vector<double>> data;
        std::stringstream ss(s);
        std::string row_str;
        while (std::getline(ss, row_str, ';')) {
            std::vector<double> row;
            std::stringstream rss(row_str);
            std::string val_str;
            while (std::getline(rss, val_str, ',')) row.push_back(std::stod(val_str));
            data.push_back(row);
        }
        int rows = (int)data.size(), cols = rows > 0 ? (int)data[0].size() : 0;
        FloatMatrix mat(rows, cols);
        for(int i=0; i<rows; ++i)
            for(int j=0; j<(int)data[i].size(); ++j) mat(i,j) = data[i][j];
        return mat;
    }

public:
    void add_constant(std::string name, TypedShape ts) { constants[name] = ts; }
    
    void add_function(std::string name, DataType ret, std::vector<DataType> args, 
                      std::function<TypedShape(const std::vector<TypedShape>&)> res, int id) {
        registry[name] = {ret, args, res, id};
    }

    // Shunting-Yard Tokenizer
    std::vector<std::string> to_rpn(const std::string& expr) {
        std::vector<std::string> output;
        std::stack<std::string> ops;
        std::regex re(R"(\.[tT]|[a-zA-Z_][a-zA-Z0-9_]*|\[[^\]]+\]|\d+\.?\d*([eE][+-]?\d+)?|\*|\+|\(|\)|,)");
        std::vector<std::string> tokens{std::sregex_token_iterator(expr.begin(), expr.end(), re), std::sregex_token_iterator()};

        for (const auto& t : tokens) {
            if (constants.count(t) || is_number(t) || t[0] == '[') output.push_back(t);
            else if (registry.count(t) || t == "(") ops.push(t);
            else if (t == ")") {
                while (!ops.empty() && ops.top() != "(") { output.push_back(ops.top()); ops.pop(); }
                if (!ops.empty()) ops.pop();
                if (!ops.empty() && registry.count(ops.top())) { output.push_back(ops.top()); ops.pop(); }
            } else if (precedence.count(t)) {
                while (!ops.empty() && precedence.count(ops.top()) && precedence[ops.top()] >= precedence[t]) {
                    output.push_back(ops.top()); ops.pop();
                }
                ops.push(t);
            }
        }
        while (!ops.empty()) { output.push_back(ops.top()); ops.pop(); }
        return output;
    }

    // Compile into Bytecode with Scalar-Aware Dimension Checks
    void compile(const std::string& expr, std::vector<Instruction>& program, 
                 std::vector<FloatMatrix>& pool, std::map<std::string, int>& symbols) {
        auto rpn = to_rpn(expr);
        std::stack<int> idxStack;
        std::stack<TypedShape> shapeStack;
        std::stack<bool> transStack; 

        for (const auto& t : rpn) {
            if (t[0] == '[') { // Matrix Literal
                FloatMatrix mat = parse_matrix_literal(t);
                idxStack.push((int)pool.size());
                shapeStack.push(TypedShape::Matrix((int)mat.rows(), (int)mat.cols()));
                pool.push_back(mat);
                transStack.push(false);
            }
            else if (is_number(t)) { // Scalar Literal
                idxStack.push((int)pool.size());
                shapeStack.push(TypedShape::Scalar());
                pool.push_back(FloatMatrix::fromScalar(std::stod(t)));
                transStack.push(false);
            }
            else if (constants.count(t)) { // Named Variable
                symbols[t] = (int)pool.size();
                idxStack.push((int)pool.size());
                shapeStack.push(constants[t]);
                pool.push_back(FloatMatrix(constants[t].rows, constants[t].cols));
                transStack.push(false);
            }
            else if (registry.count(t)) { // Function Call
                auto& sig = registry[t];
                std::vector<TypedShape> args; std::vector<int> arg_idx;
                for(size_t i=0; i<sig.argTypes.size(); ++i) {
                    args.insert(args.begin(), shapeStack.top()); shapeStack.pop();
                    arg_idx.insert(arg_idx.begin(), idxStack.top()); idxStack.pop();
                    transStack.pop(); 
                }
                TypedShape res_s = sig.resolver(args);
                int res_i = (int)pool.size();
                pool.push_back(FloatMatrix(res_s.rows, res_s.cols));
                program.push_back({OpCode::CALL_FUNC, sig.provider_id, arg_idx, -1, -1, false, false, res_i});
                idxStack.push(res_i); shapeStack.push(res_s); transStack.push(false);
            } 
            else if (t == ".T") { // Transpose View
                TypedShape s = shapeStack.top(); shapeStack.pop();
                shapeStack.push(TypedShape::Matrix(s.cols, s.rows));
                bool trans = transStack.top(); transStack.pop();
                transStack.push(!trans); 
            } 
            else if (t == "*") { // Multiplication
                bool rT = transStack.top(); transStack.pop();
                bool lT = transStack.top(); transStack.pop();
                TypedShape rs = shapeStack.top(); shapeStack.pop();
                TypedShape ls = shapeStack.top(); shapeStack.pop();
                int ri = idxStack.top(); idxStack.pop();
                int li = idxStack.top(); idxStack.pop();

                int outR, outC;
                // --- SCALAR MULTIPLICATION LOGIC ---
                if (ls.type == DataType::SCALAR) { outR = rs.rows; outC = rs.cols; }
                else if (rs.type == DataType::SCALAR) { outR = ls.rows; outC = ls.cols; }
                else { // Regular Matrix Multiply
                    if (ls.cols != rs.rows) {
                        std::cerr << "COMPILE ERROR: Mismatch " << ls.rows << "x" << ls.cols << " * " << rs.rows << "x" << rs.cols << std::endl;
                        exit(1);
                    }
                    outR = ls.rows; outC = rs.cols;
                }

                int res_i = (int)pool.size();
                pool.push_back(FloatMatrix(outR, outC));
                program.push_back({OpCode::MAT_MUL, -1, {}, li, ri, lT, rT, res_i});
                idxStack.push(res_i); shapeStack.push(TypedShape::Matrix(outR, outC)); transStack.push(false);
            }
        }
    }
};

/**
 * --- PART 4: THE EVALUATOR ---
 */

class FEAEvaluator {
    std::vector<Instruction> program;
    std::vector<FloatMatrix> pool;
    std::map<std::string, int> symbols;
    using Provider = std::function<void(const std::vector<const FloatMatrix*>&, FloatMatrix&)>;
    std::vector<Provider> providers;

public:
    FEAEvaluator(std::vector<Instruction> pr, std::vector<FloatMatrix> po, std::map<std::string, int> sy, std::vector<Provider> pv)
        : program(pr), pool(po), symbols(sy), providers(pv) {}

    void set_parameter(const std::string& name, const FloatMatrix& data) {
        if (symbols.count(name)) pool[symbols[name]] = data;
    }

    void execute() {
        for (const auto& instr : program) {
            if (instr.op == OpCode::CALL_FUNC) {
                std::vector<const FloatMatrix*> args;
                for (int i : instr.inputs) args.push_back(&pool[i]);
                providers[instr.func_id](args, pool[instr.output_idx]);
            } else {
                auto& L = pool[instr.lhs_idx]; auto& R = pool[instr.rhs_idx]; auto& O = pool[instr.output_idx];
                /* 
                // EIGEN native implementation with optimizations for scalar cases 
                // Optimized execution: detect if one operand is a 1x1 scalar
                if (L.rows() == 1 && L.cols() == 1)      O.noalias() = L(0,0) * (instr.rhs_trans ? R.transpose() : R);
                else if (R.rows() == 1 && R.cols() == 1) O.noalias() = (instr.lhs_trans ? L.transpose() : L) * R(0,0);
                else { // Standard Matrix-Matrix product
                    auto LV = instr.lhs_trans ? L.transpose() : L.eval(); // eval() simplified for clarity
                    auto RV = instr.rhs_trans ? R.transpose() : R.eval();
                    // In real implementation, use the 4-way if-else from previous version for max speed
                    if (!instr.lhs_trans && !instr.rhs_trans) O.noalias() = L * R;
                    else if (instr.lhs_trans && !instr.rhs_trans) O.noalias() = L.transpose() * R;
                    else if (!instr.lhs_trans && instr.rhs_trans) O.noalias() = L * R.transpose();
                    else O.noalias() = L.transpose() * R.transpose();
                }
                */
               // Optimized execution: detect if one operand is a 1x1 scalar
                if (L.rows() == 1 && L.cols() == 1) {
                    if (instr.rhs_trans) {O.beTranspositionOf(R);} else {O=R;}
                    O.times(L(0,0));
                } else if (R.rows() == 1 && R.cols() == 1) {
                    if (instr.lhs_trans) {O.beTranspositionOf(L);} else {O=L;}
                    O.times(R(0,0));
                }
                else { // Standard Matrix-Matrix product
                    if (!instr.lhs_trans && !instr.rhs_trans)      O.beProductOf(L, R);
                    else if (instr.lhs_trans && !instr.rhs_trans)  O.beTProductOf(L, R);
                    else if (!instr.lhs_trans && instr.rhs_trans)  O.beProductTOf(L, R);
                    else                                           O.beTProductTOf(L, R);
                }
            }
        }
    }
    const FloatMatrix& result() const { return pool.back(); }
};

} // namespace oofem
#endif // mpmevaluator_h