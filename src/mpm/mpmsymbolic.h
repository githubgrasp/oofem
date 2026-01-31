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
#ifndef mpmsymbolic_h
#define mpmsymbolic_h

/**
 * Multiphysics module 
 * Classes:
 * - ElementBase(Element) defining geometry
 * - Variable class representing unknown field (or test feld) in a weak psolution. The variable has its interpolation, type (scalar, vector), size.
    When test field, it keeps reference to its primary (unknown) variable. The history parameter dermines how many time steps to remember. 
 * - Term class represnting a term to evaluate on element. Paramaters element(geometry), variables
 * - Element - responsible for defining and performing integration (of terms), assembly of term contributions. 
 */

#include "classfactory.h"
#include "mpm.h"
#include "mpmevaluator2.h"
#include "logger.h"


namespace oofem {
/**
 * @brief Symbolic term allowing to parse and evaluate user defined expressions
 * 
 */
class SymbolicTerm : public GenericCellTerm {
    protected:
        std::string expression;
        mutable int pool_ptr=0;
        mutable std::vector<Instruction> program;
        mutable std::map<std::string, int> symbols;
        mutable std::map<int, VarData> constants;
 
        struct TestField {
            FloatMatrix values; // Nodal values (e.g., Temperature at 3 nodes)
            TestField(double start_val) { 
            values=FloatMatrix::fromIniList({{start_val}, {start_val + 10.0}, {start_val + 20.0}}); 
        }
    };

    public:
    SymbolicTerm() : GenericCellTerm() {}
    SymbolicTerm (const Variable *testField, const Variable* unknownField, const std::string& expr, MaterialMode m=MaterialMode::_Unknown)  : GenericCellTerm(testField, unknownField, m), expression(expr) {}

    void initializeFrom(const std::shared_ptr<InputRecord> &ir, EngngModel* problem) override {
        GenericCellTerm::initializeFrom(ir, problem);
        IR_GIVE_FIELD(ir, expression, "expression");

        MPMCompiler compiler;

        // Register variables (classes)
        compiler.register_function("B"); 

        try {
            compiler.compile_script(expression, program, symbols, constants, pool_ptr);
        } catch (const std::exception& e) {
            std::string msg = "SymbolicTerm: Compilation error in expression '" + expression + "': " + e.what();
            OOFEM_LOG_ERROR("%s", msg.c_str());
        }
     }
    void evaluate_lin (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override {
        try {   
            MPMEvaluator vm(pool_ptr, symbols);
            for(auto const& [idx, val] : constants) vm.init_slot(idx, val);

 
            vm.set_variable("u", 1.0);
            //vm.set_variable("v", 1.0, VarType::SCALAR);
            vm.register_functor("B", [](const auto& a, auto& out){
                FloatMatrix m = FloatMatrix::fromIniList({{1.}, {-1.}});
                std::cout << "B called\n" << m << "\n";
                out.value = m; out.type = VarSlot::Type::MATRIX;
            });
        vm.execute(program);
        if (vm.get_result().type == VarSlot::Type::MATRIX) {
            FloatMatrix m = std::get<FloatMatrix>(vm.get_result().value);
            std::cout << "Result:\n" << m << "\n";
        }

        } catch (const std::exception& e) {
            std::cerr << "VM ERROR: " << e.what() << std::endl;
        }
    }
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override{}
    void getDimensions(Element& cell) const override {}

};

}
#endif // mpmsymbolic_h