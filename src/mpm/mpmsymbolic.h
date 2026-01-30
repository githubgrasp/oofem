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
        mutable int pool_size;
        mutable std::vector<Instruction> program;
        mutable std::map<std::string, int> symbols;
        mutable std::vector<std::pair<int, FloatMatrix>> constants;
        mutable std::vector<Provider> providers;
    
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
        compiler.add_variable("u");
        compiler.add_variable("v");
        // register functors
        compiler.add_function("GetNodal", 0);
        auto Nodal_provider = [](ObjectRegistry& reg, const auto& args, auto& out) {
            if (args[0]->type != VarType::OBJECT) throw std::runtime_error("Arg 0 must be an Object");
            // 1. Recover the Symbolic IDs from the 1x1 matrices in args
            int id_field = (int)args[0]->data(0,0);    
            // 2. Resolve to C++ Pointers
            TestField* f = static_cast<TestField*>(reg[id_field]);
            if(!f) throw std::runtime_error("Dynamic resolution failed!");

            out.data = f->values; 
            out.type = VarType::MATRIX;
        };
        providers.push_back(Nodal_provider);

        //std::string expr = "(GetNodal(u) + GetNodal(v))* 1.5e2";
        
        try {
            compiler.compile(expression, program, symbols, constants, pool_size);
        } catch (const std::exception& e) {
            std::string msg = "SymbolicTerm: Compilation error in expression '" + expression + "': " + e.what();
            OOFEM_LOG_ERROR("%s", msg.c_str());
        }
     }
    void evaluate_lin (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override {
        try {   
           MPMEvaluator vm(program, pool_size, constants, providers);
            TestField fu(10.0); // Example nodal values
            TestField fv(100.0);

            vm.bind_object(this->symbols["u"], &fu);
            vm.bind_object(this->symbols["v"], &fv);
            
            vm.execute();
            std::cout << "Success! Calculated Scalar Result:\n" << vm.result() << std::endl;
            answer= vm.result();

        } catch (const std::exception& e) {
            std::cerr << "VM ERROR: " << e.what() << std::endl;
        }

    }
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override{}
    void getDimensions(Element& cell) const override {}

};

}
#endif // mpmsymbolic_h