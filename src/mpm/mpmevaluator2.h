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
 * --- DATA LAYER ---
 */

enum class OpCode { CALL_FUNC, MAT_MUL, MAT_ADD, MAT_SUB, MAT_NEG, MAT_DOT };

struct Instruction {
    OpCode op;
    int func_id;             
    std::vector<int> inputs; 
    int lhs_idx = -1;        
    int rhs_idx = -1;        
    bool lhs_trans = false; 
    bool rhs_trans = false;
    int output_idx;      
};

/**
 * --- COMPILER LAYER ---
 */
class MPMCompiler {
    // PRECEDENCE TABLE
    // 6: Unary minus (-u)
    // 5: Methods (.dot, .T) - High precedence ensures scalar reduction happens first
    // 4: Multiplication (*)
    // 1: Addition/Subtraction (+, -)
    std::map<std::string, int> precedence = {
        {"unary-", 6}, {".T", 5}, {"*", 4}, {"+", 1}, {"-", 1}
    };
    std::map<std::string, int> function_registry;
    std::vector<std::string> variables;

    FloatMatrix parse_literal(std::string s) {
        if (s[0] != '[') return FloatMatrix::fromScalar(std::stod(s));
        
        s = std::regex_replace(s, std::regex(R"([\[\]])"), "");
        std::replace(s.begin(), s.end(), ',', ' ');

        std::vector<std::vector<double>> data;
        std::stringstream ss(s);
        std::string row_str;

        while (std::getline(ss, row_str, ';')) {
            if (row_str.find_first_not_of(" \t\n\r") == std::string::npos) continue;
            std::vector<double> row;
            std::stringstream rss(row_str);
            double val;
            while (rss >> val) row.push_back(val);
            if (!row.empty()) data.push_back(row);
        }

        int rows = (int)data.size(), cols = (rows > 0) ? (int)data[0].size() : 0;
        FloatMatrix mat(rows, cols);
        for(int i=0; i<rows; ++i)
            for(int j=0; j<cols; ++j) mat(i,j) = data[i][j];
        return mat;
    }

public:
    void add_function(std::string name, int id) { function_registry[name] = id; }
    void add_variable(std::string name) { variables.push_back(name); }

    std::vector<std::string> to_rpn(const std::string& expr) {
        std::vector<std::string> output;
        std::stack<std::string> ops;
        std::stack<int> arg_counts;

        // REGEX: Enhanced to capture scientific notation with optional signs (e.g. 1e-10)
        std::regex re(R"(\.dot|\.[tT]|[a-zA-Z_][a-zA-Z0-9_]*|\[[^\]]+\]|\d*\.?\d+(?:[eE][-+]?\d+)?|\*|\+|\-|\(|\)|,)");
        auto it = std::sregex_token_iterator(expr.begin(), expr.end(), re);
        auto end = std::sregex_token_iterator();
        std::vector<std::string> tokens(it, end);

        std::string last_token = "(";
        for (const auto& t : tokens) {
            // Check for Leaf nodes: Numbers, Literals, or registered variables
            if (isdigit(t[0]) || (t[0] == '.' && t.size() > 1 && isdigit(t[1])) || t[0] == '[' || 
                std::find(variables.begin(), variables.end(), t) != variables.end()) {
                output.push_back(t);
                if (!arg_counts.empty()) arg_counts.top()++;
            } 
            else if (function_registry.count(t)) {
                ops.push(t); arg_counts.push(0);
            } 
            else if (t == "(") {
                ops.push(t);
            } 
            else if (t == ",") {
                while (!ops.empty() && ops.top() != "(") { output.push_back(ops.top()); ops.pop(); }
            } 
            else if (t == ")") {
                while (!ops.empty() && ops.top() != "(") { output.push_back(ops.top()); ops.pop(); }
                if(!ops.empty()) ops.pop(); // Pop '('
                if (!ops.empty() && function_registry.count(ops.top())) {
                    output.push_back(ops.top() + "#" + std::to_string(arg_counts.top()));
                    ops.pop(); arg_counts.pop();
                }
            } 
            else if (precedence.count(t)) {
                std::string op = t;
                // Unary minus detection: if '-' follows an operator or bracket, it's unary
                if (op == "-" && (last_token == "(" || last_token == "," || precedence.count(last_token))) op = "unary-";
                while (!ops.empty() && ops.top() != "(" && precedence.count(ops.top()) && precedence[ops.top()] >= precedence[op]) {
                    output.push_back(ops.top()); ops.pop();
                }
                ops.push(op);
            }
            last_token = t;
        }
        while (!ops.empty()) { output.push_back(ops.top()); ops.pop(); }
        return output;
    }

    void compile(const std::string& expr, std::vector<Instruction>& program, 
                 std::map<std::string, int>& symbols, std::vector<std::pair<int, FloatMatrix>>& constants, int& pool_size) {
        auto rpn = to_rpn(expr);
        int next_idx = 0;
        std::stack<int> idxStack;
        std::stack<bool> transStack;

        for (const auto& t : rpn) {
            if (t.find('#') != std::string::npos) {
                size_t sep = t.find('#');
                int num_args = std::stoi(t.substr(sep + 1));
                std::vector<int> args(num_args);
                for (int i = num_args - 1; i >= 0; --i) { args[i] = idxStack.top(); idxStack.pop(); transStack.pop(); }
                int res = next_idx++;
                program.push_back({OpCode::CALL_FUNC, function_registry[t.substr(0, sep)], args, -1, -1, false, false, res});
                idxStack.push(res); transStack.push(false);
            }  
            else if (t == ".T") {
                bool tr = transStack.top(); transStack.pop();
                transStack.push(!tr);
            } 
            else if (precedence.count(t) || t == "unary-") {
                bool rT = transStack.top(); transStack.pop();
                int ri = idxStack.top(); idxStack.pop();
                int li = -1, res = next_idx++;  
                bool lT = false;
                if (t != "unary-") { lT = transStack.top(); transStack.pop(); li = idxStack.top(); idxStack.pop(); }
                OpCode op = (t == "*") ? OpCode::MAT_MUL : (t == "+") ? OpCode::MAT_ADD : (t == "-") ? OpCode::MAT_SUB : OpCode::MAT_NEG;
                program.push_back({op, -1, {}, li, ri, lT, rT, res});
                idxStack.push(res); transStack.push(false);
            } 
            else {
                if (symbols.find(t) == symbols.end()) {
                    symbols[t] = next_idx++;
                    if (isdigit(t[0]) || t[0] == '.' || t[0] == '[') 
                        constants.push_back({symbols[t], parse_literal(t)});
                }
                idxStack.push(symbols[t]); transStack.push(false);
            }
        }
        pool_size = next_idx;
    }
};

/**
 * --- EVALUATOR LAYER ---
 */
class MPMEvaluator {
    std::vector<Instruction> program;
    std::vector<FloatMatrix> pool;
    std::map<std::string, void*> class_registry;
    using Provider = std::function<void(std::map<std::string, void*>&, const std::vector<const FloatMatrix*>&, FloatMatrix&)>;
    std::vector<Provider> providers;

public:
    MPMEvaluator(std::vector<Instruction> pr, int sz, const std::vector<std::pair<int, FloatMatrix>>& constants, std::vector<Provider> pvs)
        : program(pr), providers(pvs) {
        pool.resize(sz);
        for (auto& c : constants) pool[c.first] = c.second;
    }

    void bind_class(const std::string& name, void* ptr) { class_registry[name] = ptr; }

    void execute() {
        for (size_t i = 0; i < program.size(); ++i) {
            const auto& instr = program[i];
            auto& O = pool[instr.output_idx];
            const auto& R = pool[instr.rhs_idx];
            const auto& L = pool[instr.lhs_idx];
 
            if (instr.op == OpCode::CALL_FUNC) {
                std::vector<const FloatMatrix*> args;
                for (int idx : instr.inputs) args.push_back(&pool[idx]);
                providers[instr.func_id](class_registry, args, O);
                continue;
            } else if (instr.op == OpCode::MAT_NEG) { 
                if (instr.rhs_trans) O.beTranspositionOf(R);
                else O = R;
                continue; 
            } else if (instr.op == OpCode::MAT_MUL) {
                // Implicit scalar broadcasting
                // Optimized execution: detect if one operand is a 1x1 scalar
                if (L.rows() == 1 && L.cols() == 1) {
                    if (instr.rhs_trans) {O.beTranspositionOf(R);} else {O=R;}
                    O.times(L(0,0));
                } else if (R.rows() == 1 && R.cols() == 1) {
                    if (instr.lhs_trans) {O.beTranspositionOf(L);} else {O=L;}
                    O.times(R(0,0));
                }
                else { // Standard Matrix-Matrix product
                    // check dimensions
                    int lcols = (instr.lhs_trans) ? L.giveNumberOfRows() : L.giveNumberOfColumns();
                    int rrows = (instr.rhs_trans) ? R.giveNumberOfColumns() : R.giveNumberOfRows();
                    if (lcols != rrows) {
                        throw std::runtime_error("Step " + std::to_string(i) + " Multiplication Mismatch: LHS Cols(" + std::to_string(lcols) + ") != RHS Rows(" + std::to_string(rrows) + ")");
                    }   
                    if (!instr.lhs_trans && !instr.rhs_trans)      O.beProductOf(L, R);
                    else if (instr.lhs_trans && !instr.rhs_trans)  O.beTProductOf(L, R);
                    else if (!instr.lhs_trans && instr.rhs_trans)  O.beProductTOf(L, R);
                    else                                           O.beTProductTOf(L, R);
                }
            } else if (instr.op == OpCode::MAT_ADD) {
                //O.noalias() = LV + RV; break;
                // check dimensions
                int lcols = (instr.lhs_trans) ? L.giveNumberOfRows() : L.giveNumberOfColumns();
                int lrows = (instr.lhs_trans) ? L.giveNumberOfColumns() : L.giveNumberOfRows();
                int rcols = (instr.rhs_trans) ? R.giveNumberOfRows() : R.giveNumberOfColumns();
                int rrows = (instr.rhs_trans) ? R.giveNumberOfColumns() : R.giveNumberOfRows();
                if (lrows != rrows || lcols != rcols) {
                    throw std::runtime_error("Step " + std::to_string(i) + " ADD Mismatch: LHS(" + std::to_string(lrows) + ", " + std::to_string(lcols) + ") != RHS (" + std::to_string(rrows) + ", " + std::to_string(rcols) + ")");
                }   
                if (instr.lhs_trans) {
                    O.beTranspositionOf(L);
                } else {
                    O = L;
                }
                if (instr.rhs_trans) {
                    FloatMatrix _RT; _RT.beTranspositionOf(R); O.add(_RT);
                } else {
                    O.add(R);
                }
                
            } else if (instr.op == OpCode::MAT_SUB) {
                //O.noalias() = LV - RV; break;
                int lcols = (instr.lhs_trans) ? L.giveNumberOfRows() : L.giveNumberOfColumns();
                int lrows = (instr.lhs_trans) ? L.giveNumberOfColumns() : L.giveNumberOfRows();
                int rcols = (instr.rhs_trans) ? R.giveNumberOfRows() : R.giveNumberOfColumns();
                int rrows = (instr.rhs_trans) ? R.giveNumberOfColumns() : R.giveNumberOfRows();
                if (lrows != rrows || lcols != rcols) {
                    throw std::runtime_error("Step " + std::to_string(i) + " SUBTRACT Mismatch: LHS(" + std::to_string(lrows) + ", " + std::to_string(lcols) + ") != RHS (" + std::to_string(rrows) + ", " + std::to_string(rcols) + ")");
                }   
                if (instr.lhs_trans) {
                    O.beTranspositionOf(L);
                } else {
                    O = L;
                }
                if (instr.rhs_trans) {
                    FloatMatrix _RT; _RT.beTranspositionOf(R); O.subtract(_RT);
                } else {
                    O.subtract(R);
                }
            }
        }// end for
    }
    const FloatMatrix& result() const { return pool.back(); }
};



} // namespace oofem
#endif // mpmevaluator_h