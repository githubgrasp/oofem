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
#include <variant>
#include <map>
#include <functional>
#include <stack>
#include <regex>
#include <sstream>
#include "floatmatrix.h"

using namespace std;

namespace oofem {
/*
The FEA Virtual Machine is a specialized execution engine designed to bridge the gap between high-level mathematical notation and high-performance C++ execution. 
Below is a summary of its core features.

1. Advanced Matrix Algebra Syntax

    The system prioritizes engineering-standard notation to ensure that expressions remain readable and maintainable.

    * Modern Transpose (`.T`): Uses the Python/NumPy convention for transposes, preventing conflicts with potential exponentiation (`^`) operators.
    * Postfix Slicing: Supports MATLAB-style slicing (e.g., `K[0, :]` or `M[:, 1]`), allowing for efficient extraction of specific degrees of freedom or sub-matrices.
    * Shape Promotion: Automatically converts  matrix results into scalars, enabling seamless use in logical and conditional expressions.

2. High-Performance Runtime Architecture

    Designed for the iterative nature of Finite Element Analysis, the runtime minimizes overhead during evaluation.

    * Static Memory Pooling: All intermediate calculation results are assigned fixed slots in a pre-allocated memory pool during compilation. This eliminates expensive heap allocations during the execution loop.
    * Broadcasting Logic: Automatically handles operations between mismatched ranks (e.g., adding a scalar to every element of a matrix) without requiring manual loops from the user.

3. Comprehensive Operator Suite

    The VM supports a full range of operators necessary for non-linear and conditional physics:

    * Arithmetic: Unary plus/minus, addition, subtraction, and matrix multiplication.
    * Relational: Comparison operators (`==`, `>`, `<`) used for identifying material state changes or boundary conditions.
    * Logical: Boolean operators (`&&`, `||`) to construct complex multi-condition branching.



4. Extensible Functor System

    The engine is not limited to its built-in operators.

    * Host-Registered Functors: Developers can register custom C++ functions (like `if()`, `sin()`, or specific element stiffness kernels) that the VM can call directly using arguments from the memory pool.
    * Multivariate Support: The compiler handles functions with a dynamic number of arguments, enabling complex branching and multi-parameter material laws.

5. Automated Compilation Pipeline

    The transformation from string to machine code is fully handled by a robust two-stage parser:

    * Shunting-Yard Parser: Resolves operator precedence and associativity to ensure mathematical correctness (e.g., transposes bind tighter than multiplications).
    * Bytecode Generator: Produces a linear stream of `Instruction` objects, making the execution phase a simple, high-speed dispatch loop.

*/


/** 
 * --- FEA VM CORE TYPES ---
 * Dynamic typing allows matrices to change shape during runtime, 
 * which is essential for element-specific calculations.
 */

/**
 * @enum VarType
 * @brief Supported data types within the VM pool.
 */
enum class VarType { SCALAR, MATRIX };
using VarData = std::variant<double, FloatMatrix>;

/**
 * @struct VarSlot
 * @brief Atomic memory unit in the VM's pre-allocated pool.
 */
struct VarSlot {
    VarType type = VarType::SCALAR;
    VarData value = 0.0;
};

using Functor = std::function<void(const std::vector<const VarSlot*>&, VarSlot&)>;

/**
 * @enum OpCode
 * @brief Operation codes for the VM instruction set.
 */
enum class OpCode { 
    CALL_FUNC, MAT_SLICE, MAT_ADD, MAT_SUB, MAT_MUL, 
    MAT_NEG, MAT_POS, MAT_TRANS, LOGIC_OP 
};

struct Instruction {
    OpCode op;
    std::string func_name; 
    std::vector<int> inputs; 
    int lhs_idx = -1, rhs_idx = -1, output_idx = -1;
};

/**
 * @class MPMCompiler
 * @brief Parses infix strings and emits optimized Bytecode.
 * * The compiler utilizes a Shunting-Yard algorithm to resolve operator precedence.
 * It maps variables and constants to fixed memory slots to eliminate runtime 
 * allocations during the execution phase.
 * * 
 */
class MPMCompiler {
    // .T is a postfix operator with the highest precedence.
    std::map<std::string, int> precedence = {
        {".T", 10}, {"unary-", 9}, {"unary+", 9}, {"*", 8}, {"+", 7}, {"-", 7},
        {"==", 6}, {">", 6}, {"<", 6}, {"&&", 5}, {"||", 4}
    };
    std::map<std::string, bool> func_registry;

    bool is_operand(const std::string& t) {
        if (t.empty()) return false;
        if (t == ")" || t == ".T" || t[0] == ']' || isalnum(t[0])) return true;
        return false;
    }

    public:
    /**
     * @brief Registers a function name for the parser to recognize.
     * @param name The string name of the function (e.g., "K", "if", "sin").
     */
    void register_function(std::string name) { func_registry[name] = true; }
    /**
     * @brief Compiles an infix math expression into VM instructions.
     * * @param expr The mathematical string (e.g., "u.T * K * u").
     * @param program [out] Vector to be filled with generated instructions.
     * @param symbols [in/out] Mapping of variable names to pool indices.
     * @param constants [out] Mapping of pool indices to literal values found in string.
     * @param pool_ptr [in/out] Current tail of the memory pool.
     */
    void compile(std::string expr, std::vector<Instruction>& program, 
                 std::map<std::string, int>& symbols, std::map<int, VarData>& constants, int& pool_ptr) {

        // Regex includes '.T' as a specific multi-character token
        // Updated regex to capture potential matrix literals and slices
        std::regex re(R"(\[\[.*?\]\]|\[.*?\]|[a-zA-Z_]\w*|\d*\.?\d+|&&|\|\||==|\.T|\+|\-|\*|\^|>|<|\(|\)|,)");
        auto it = std::sregex_token_iterator(expr.begin(), expr.end(), re);
        std::vector<std::string> tokens(it, std::sregex_token_iterator());

        std::vector<std::string> rpn;
        std::stack<std::string> ops;
        std::stack<int> arg_count_stack; 
        std::string last = "(";

        for (const std::string& t : tokens) {
            if (t[0] == '[') {
                if (is_operand(last)) parse_slice(t, rpn, constants, pool_ptr);
                else parse_literal(t, rpn, constants, pool_ptr);
            }
            else if (func_registry.count(t)) { ops.push(t); arg_count_stack.push(0); }
            else if (t == "(") { ops.push(t); }
            else if (t == ")") {
                while (!ops.empty() && ops.top() != "(") { rpn.push_back(ops.top()); ops.pop(); }
                ops.pop(); 
                if (!ops.empty() && func_registry.count(ops.top())) {
                    rpn.push_back(ops.top() + "#" + std::to_string(arg_count_stack.top()));
                    arg_count_stack.pop(); ops.pop();
                }
            }
            else if (t == ",") {
                while (!ops.empty() && ops.top() != "(") { rpn.push_back(ops.top()); ops.pop(); }
                if (!arg_count_stack.empty()) arg_count_stack.top()++;
            }
            else if (precedence.count(t)) {
                std::string op = t;
                // Distinguish between binary and unary +/-
                if (t == "-" && (last == "(" || last == ",")) op = "unary-";
                if (t == "+" && (last == "(" || last == ",")) op = "unary+";
                
                while (!ops.empty() && ops.top() != "(" && precedence.count(ops.top()) && precedence[ops.top()] >= precedence[op]) {
                    rpn.push_back(ops.top()); ops.pop();
                }
                ops.push(op);
            }
            else { 
                rpn.push_back(t); 
                if (!arg_count_stack.empty() && arg_count_stack.top() == 0) arg_count_stack.top() = 1;
            }
            last = t;
        }
        while (!ops.empty()) { rpn.push_back(ops.top()); ops.pop(); }
        emit_bytecode(rpn, program, symbols, constants, pool_ptr);
    }

private:
    void parse_literal(std::string t, std::vector<std::string>& rpn, std::map<int, VarData>& constants, int& pool_ptr) {
        std::string content = t.substr(1, t.size() - 2); 
        std::vector<std::vector<double>> rows;
        std::regex row_re(R"(\[(.*?)\])");
        auto row_it = std::sregex_iterator(content.begin(), content.end(), row_re);
        
        for (std::sregex_iterator i = row_it; i != std::sregex_iterator(); ++i) {
            std::stringstream ss((*i)[1].str());
            std::string val; std::vector<double> r_vals;
            while (std::getline(ss, val, ',')) r_vals.push_back(std::stod(val));
            rows.push_back(r_vals);
        }

        FloatMatrix m(rows.size(), rows[0].size());
        for(std::size_t i=0; i<rows.size(); ++i)
            for(std::size_t j=0; j<rows[0].size(); ++j) m(i,j) = rows[i][j];

        int slot = pool_ptr++; constants[slot] = m;
        rpn.push_back("__CONST_MAT__" + std::to_string(slot));
    }
    void parse_slice(std::string t, std::vector<std::string>& rpn, std::map<int, VarData>& constants, int& pool_ptr) {
        std::string content = t.substr(1, t.size() - 2);
        std::stringstream ss(content);
        std::string seg; int count = 0;
        while (std::getline(ss, seg, ',')) {
            seg.erase(0, seg.find_first_not_of(" "));
            if (seg == ":") {
                int slot = pool_ptr++; constants[slot] = -1.0;
                rpn.push_back("__CONST__" + std::to_string(slot));
            } else { rpn.push_back(seg); }
            count++;
        }
        rpn.push_back("SLICE#" + std::to_string(count));
    }

    void emit_bytecode(std::vector<std::string>& rpn, std::vector<Instruction>& program, 
                       std::map<std::string, int>& symbols, std::map<int, VarData>& constants, int& pool_ptr) {
        std::stack<int> eval_stack;
        for (auto& t : rpn) {
            size_t h_pos = t.find('#');
            if (h_pos != std::string::npos) {
                std::string name = t.substr(0, h_pos);
                int n = std::stoi(t.substr(h_pos + 1));
                std::vector<int> args(n);
                for (int i = n - 1; i >= 0; --i) { args[i] = eval_stack.top(); eval_stack.pop(); }
                int out = pool_ptr++;
                if (name == "SLICE") {
                    int target = eval_stack.top(); eval_stack.pop();
                    std::vector<int> in = {target}; in.insert(in.end(), args.begin(), args.end());
                    program.push_back({OpCode::MAT_SLICE, "", in, -1, -1, out});
                } else {
                    program.push_back({OpCode::CALL_FUNC, name, args, -1, -1, out});
                }
                eval_stack.push(out);
            } else if (precedence.count(t)) {
                int r = eval_stack.top(); eval_stack.pop();
                int out = pool_ptr++;
                if (t == ".T") program.push_back({OpCode::MAT_TRANS, "", {}, -1, r, out});
                else if (t == "unary-") program.push_back({OpCode::MAT_NEG, "", {}, -1, r, out});
                else if (t == "unary+") program.push_back({OpCode::MAT_POS, "", {}, -1, r, out});
                else {
                    int l = eval_stack.top(); eval_stack.pop();
                    OpCode op = (t == "*") ? OpCode::MAT_MUL : (t == "+") ? OpCode::MAT_ADD : (t == "-") ? OpCode::MAT_SUB : OpCode::LOGIC_OP;
                    program.push_back({op, t, {}, l, r, out});
                }
                eval_stack.push(out);
            } else {
                int slot;
                if (t.find("__CONST_MAT__") == 0) slot = std::stoi(t.substr(13));
                else if (t.find("__CONST__") == 0) slot = std::stoi(t.substr(9));
                else if (isdigit(t[0]) || (t.size() > 1 && t[0] == '.')) { slot = pool_ptr++; constants[slot] = std::stod(t); }
                else {
                    if (symbols.find(t) == symbols.end()) symbols[t] = pool_ptr++;
                    slot = symbols[t];
                }
                eval_stack.push(slot);
            }
        }
        // FINAL MOVE: Ensure single-atom expressions populate the final pool slot
        if (!eval_stack.empty()) {
            int final_val = eval_stack.top(); eval_stack.pop();
            int out = pool_ptr++;
            program.push_back({OpCode::MAT_POS, "", {}, -1, final_val, out});
            eval_stack.push(out);
        }
    }
};


/**
 * @class MPMEvaluator
 * @brief The execution engine that processes Bytecode instructions.
 * *The Evaluator maintains the VarSlot pool and handles the mathematical kernels 
 * via OOFEM routines. It supports "broadcasting," where a scalar can be added/subtracted 
 * from a matrix element-wise automatically.
 * *
 */
class MPMEvaluator {
    std::vector<VarSlot> pool;
    std::map<std::string, int> symbol_table;
    std::map<std::string, Functor> functors;
    std::vector<bool> initialized; // Tracks if a slot has been written to

public:
    /**
     * @brief Constructs an evaluator with a fixed memory footprint.
     * @param sz Total number of slots required for variables, constants, and intermediates.
     * @param syms The symbol table generated during compilation.
     */
    MPMEvaluator(int sz, std::map<std::string, int> syms) : pool(sz), symbol_table(syms), initialized(sz, false) {}
    /**
     * @brief Injects a C++ function into the VM environment.
     * @param name Name of the function as it appears in the expression string.
     * @param f Lambda or function pointer following the Functor signature.
     */
    void register_functor(std::string name, Functor f) { functors[name] = f; }
    void set_variable(std::string name, VarData val, VarType type) { int idx = symbol_table.at(name); pool[idx] = {type, val}; initialized[idx] = true; /* Mark as safe to read */ }
    void init_constant(int idx, VarData d, VarType t) { pool[idx] = {t, d}; initialized[idx] = true; }
    /**
     * @brief Runs the compiled instruction set.
     * @param program Vector of instructions to execute.
     * @throw std::runtime_error if dimension mismatches occur during matrix ops.
     */
    void execute(const std::vector<Instruction>& program) {
        for (const auto& instr : program) {
            // Check inputs for initialization
            if (instr.lhs_idx != -1 && !initialized[instr.lhs_idx]) 
                throw std::runtime_error("Runtime Error: Variable at slot " + std::to_string(instr.lhs_idx) + " is uninitialized.");
            if (instr.rhs_idx != -1 && !initialized[instr.rhs_idx])
                throw std::runtime_error("Runtime Error: Variable at slot " + std::to_string(instr.rhs_idx) + " is uninitialized.");
            for(int in : instr.inputs)
                if(!initialized[in]) throw std::runtime_error("Runtime Error: Function argument at slot " + std::to_string(in) + " is uninitialized.");

            VarSlot& O = pool[instr.output_idx];
            switch (instr.op) {
                case OpCode::MAT_ADD: 
                case OpCode::MAT_SUB: {
                    const VarSlot &L = pool[instr.lhs_idx], &R = pool[instr.rhs_idx];
                    bool is_add = (instr.op == OpCode::MAT_ADD);
                    if (L.type == VarType::SCALAR && R.type == VarType::SCALAR) {
                        O.value = is_add ? std::get<double>(L.value) + std::get<double>(R.value) : std::get<double>(L.value) - std::get<double>(R.value);
                        O.type = VarType::SCALAR;
                    } else if (L.type == VarType::MATRIX && R.type == VarType::MATRIX) {
                        const auto &m1 = std::get<FloatMatrix>(L.value), &m2 = std::get<FloatMatrix>(R.value);
                        if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) throw std::runtime_error("Addition Dimension Mismatch");
                        O.value = is_add ? FloatMatrix(m1 + m2) : FloatMatrix(m1 - m2);
                        O.type = VarType::MATRIX;
                    } else { // Broadcasting
                        throw std::runtime_error("Matrix scalar addition/substraction not supported");
                        //double s = (L.type == VarType::SCALAR) ? std::get<double>(L.value) : std::get<double>(R.value);
                        //const auto& m = (L.type == VarType::MATRIX) ? std::get<FloatMatrix>(L.value) : std::get<FloatMatrix>(R.value);
                        //if (!is_add && L.type == VarType::SCALAR) O.value = FloatMatrix(s - m.array());
                        //else O.value = is_add ? FloatMatrix(m + s) : FloatMatrix(m - s);
                        //O.type = VarType::MATRIX;
                    }
                    break;
                }
                case OpCode::MAT_MUL: {
                    const VarSlot &L = pool[instr.lhs_idx], &R = pool[instr.rhs_idx];
                    if (L.type == VarType::SCALAR && R.type == VarType::MATRIX) {
                        O.value = FloatMatrix(std::get<double>(L.value) * std::get<FloatMatrix>(R.value));
                        O.type = VarType::MATRIX;
                    } else if (L.type == VarType::MATRIX && R.type == VarType::MATRIX) {
                        const auto &m1 = std::get<FloatMatrix>(L.value), &m2 = std::get<FloatMatrix>(R.value);
                        if (m1.cols() != m2.rows()) throw std::runtime_error("Matrix Multiplication size mismatch");
                        auto res = m1 * m2;
                        if (res.rows() == 1 && res.cols() == 1) { O.value = res(0,0); O.type = VarType::SCALAR; }
                        else { O.value = res; O.type = VarType::MATRIX; }
                    } else { O.value = std::get<double>(L.value) * std::get<double>(R.value); O.type = VarType::SCALAR; }
                    break;
                }
                case OpCode::MAT_TRANS: {
                    const auto& src = pool[instr.rhs_idx];
                    if (src.type == VarType::SCALAR) O = src;
                    else { O.value=FloatMatrix::fromMatrix(std::get<FloatMatrix>(src.value), true); O.type = VarType::MATRIX; }
                    break;
                }
                case OpCode::MAT_SLICE: {
                    const auto& m = std::get<FloatMatrix>(pool[instr.inputs[0]].value);
                    int r = (int)std::get<double>(pool[instr.inputs[1]].value);
                    int c = (int)std::get<double>(pool[instr.inputs[2]].value);
                    if (r == -1) { O.value = FloatMatrix(m.col(c)); O.type = VarType::MATRIX; }
                    else if (c == -1) { O.value = FloatMatrix(m.row(r)); O.type = VarType::MATRIX; }
                    else { O.value = m(r, c); O.type = VarType::SCALAR; }
                    break;
                }
                case OpCode::MAT_NEG: {
                    const auto& src = pool[instr.rhs_idx];
                    if (src.type == VarType::SCALAR) { O.value = -std::get<double>(src.value); O.type = VarType::SCALAR; }
                    else { O.value = FloatMatrix(std::get<FloatMatrix>(src.value)); std::get<FloatMatrix>(O.value).negated(); O.type = VarType::MATRIX; }
                    break;
                }
                case OpCode::LOGIC_OP: {
                    double l = std::get<double>(pool[instr.lhs_idx].value), r = std::get<double>(pool[instr.rhs_idx].value);
                    bool res = (instr.func_name == "==") ? std::abs(l-r) < 1e-9 : (instr.func_name == "&&") ? (l!=0 && r!=0) : (instr.func_name == ">") ? l > r : l < r;
                    O.value = res ? 1.0 : 0.0; O.type = VarType::SCALAR; break;
                }
                case OpCode::CALL_FUNC: {
                    std::vector<const VarSlot*> args; for (int idx : instr.inputs) args.push_back(&pool[idx]);
                    functors.at(instr.func_name)(args, O); break;
                }
                case OpCode::MAT_POS: { O = pool[instr.rhs_idx]; break; }
            }
            initialized[instr.output_idx] = true; // Mark result as initialized for next steps
        }
    }
    /**
     * @brief Retrieves the final value from the last executed instruction.
     * @return The VarSlot containing the final result.
     */
    VarSlot get_result() { return pool.back(); }
};



} // namespace oofem
#endif // mpmevaluator_h