#include "mpmevaluator2.h"
#include "floatmatrix.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <cmath>
using namespace oofem;

TEST_CASE( "Testing mpm evaluator", "[MPMCompiler, MPMEvaluator]" ) {
    // This setup will be done 4 times in total, once for each section
    MPMCompiler compiler;
    std::vector<Instruction> prog;
    std::map<std::string, int> syms;
    std::map<int, VarData> consts;
    int ptr = 0;

    SECTION( "Simple constant expression 1" ) {
        compiler.compile_script("5 + 3 * 2", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 11.0) < 1e-9 );
    }
    SECTION( "Simple constant expression 2" ) {
        compiler.compile_script("10 - 2*3 ", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 4.0) < 1e-9 );
    }
    SECTION( "Simple constant expression 3" ) {
        compiler.compile_script("4 * (2 + 3)", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 20.0) < 1e-9 );
    }
    SECTION( "Simple constant expression 4" ) {
        compiler.compile_script("-2*(2-6)", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 8.0) < 1e-9 );
        REQUIRE( val > 0.0 );
    }
    SECTION( "Inline matrix expression 1" ) {
        compiler.compile_script("[[1,2],[3,4]] + [[5,6],[7,8]]", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 2 );
        REQUIRE( std::abs(mat(0,0) - 6.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,1) - 12.0) < 1e-9 );
    }
    SECTION( "Inline matrix expression 2" ) {
        compiler.compile_script("[[2,4,6],[8,10,12]] * 0.5", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 3 );
        REQUIRE( std::abs(mat(0,2) - 3.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,0) - 4.0) < 1e-9 );
    }
    SECTION("Inline matrix expression 3") {
        compiler.compile_script("[[1,2],[3,4]] * [[5,6],[7,8]]", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 2 );
        REQUIRE( std::abs(mat(0,0) - 19.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,1) - 50.0) < 1e-9 );
    }
    SECTION("Matrix slicing 1") {
        compiler.compile_script("A = [[10,20,30],[40,50,60]]; A[1,1] + A[1,2]", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 110.0) < 1e-9 );
    }
    SECTION("Matrix slicing 2") {
        compiler.compile_script("A = [[1,2,3],[4,5,6]]; A[:,1]", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 1 );
        REQUIRE( std::abs(mat(0,0) - 2.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,0) - 5.0) < 1e-9 );
    }
    SECTION("Matrix transpose") {
        compiler.compile_script("A = [[1,2,3],[4,5,6]]; A.T", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 3 );
        REQUIRE( mat.cols() == 2 );
        REQUIRE( std::abs(mat(0,1) - 4.0) < 1e-9 );
        REQUIRE( std::abs(mat(2,0) - 3.0) < 1e-9 );
    }
    SECTION("Mapped variable expression") {
        compiler.compile_script("C = A * B", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        // Setup variables A and B
        FloatMatrix A(2,2); A(0,0)=1.0; A(0,1)=2.0; A(1,0)=3.0; A(1,1)=4.0; 
        FloatMatrix B(2,2); B(0,0)=5.0; B(0,1)=6.0; B(1,0)=7.0; B(1,1)=8.0;
        vm.set_variable("A", A);
        vm.set_variable("B", B);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 2 );
        REQUIRE( std::abs(mat(0,0) - 19.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,1) - 50.0) < 1e-9 );
    }
    SECTION("Error on uninitialized variable") {
        compiler.compile_script("C = A + B", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        // Setup only variable A
        FloatMatrix A(2,2); A(0,0)=1.0; A(0,1)=2.0; A(1,0)=3.0; A(1,1)=4.0; 
        vm.set_variable("A", A);
        REQUIRE_THROWS_AS( vm.execute(prog), std::runtime_error );
    }   
    SECTION("Error on dimension mismatch") {
        compiler.compile_script("C = A * B", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        // Setup variables A and B with incompatible sizes
        FloatMatrix A(2,3); A(0,0)=1.0; A(0,1)=2.0; A(0,2)=3.0; A(1,0)=4.0; A(1,1)=5.0; A(1,2)=6.0; 
        FloatMatrix B(2,2); B(0,0)=5.0; B(0,1)=6.0; B(1,0)=7.0; B(1,1)=8.0;
        vm.set_variable("A", A);
        vm.set_variable("B", B);
        REQUIRE_THROWS_AS( vm.execute(prog), std::runtime_error );
    }
    SECTION("Function call") {
        // Register a simple functor that adds two matrices
        auto add_matrices = [](const std::vector<const VarSlot*>& inputs, VarSlot& output) {
            const FloatMatrix& A = std::get<FloatMatrix>(inputs[0]->value);
            const FloatMatrix& B = std::get<FloatMatrix>(inputs[1]->value);
            if (A.rows() != B.rows() || A.cols() != B.cols()) throw std::runtime_error("Dimension mismatch in add_matrices");
            FloatMatrix C(A.rows(), A.cols());
            for (int r = 0; r < A.rows(); ++r)
                for (int c = 0; c < A.cols(); ++c)
                    C(r,c) = A(r,c) + B(r,c);
            output.value = FloatMatrix(C);
            output.type = VarSlot::Type::MATRIX;
        };
        compiler.register_function("add"); 
        compiler.compile_script("add(A, B)", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.register_functor("add", add_matrices);
        // Setup variables A and B
        FloatMatrix A(2,2); A(0,0)=1.0; A(0,1)=2.0; A(1,0)=3.0; A(1,1)=4.0; 
        FloatMatrix B(2,2); B(0,0)=5.0; B(0,1)=6.0; B(1,0)=7.0; B(1,1)=8.0;
        vm.set_variable("A", A);
        vm.set_variable("B", B);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        FloatMatrix mat = std::get<FloatMatrix>(res.value);
        REQUIRE( mat.rows() == 2 );
        REQUIRE( mat.cols() == 2 );
        REQUIRE( std::abs(mat(0,0) - 6.0) < 1e-9 );
        REQUIRE( std::abs(mat(1,1) - 12.0) < 1e-9 );
    }
    SECTION ("Logical operations") {
        compiler.compile_script("res = (a > b) && (c < d)", prog, syms, consts, ptr);
        MPMEvaluator vm(ptr, syms);
        for(auto const& [idx, val] : consts) vm.init_slot(idx, val);
        vm.set_variable("a", 5.0);
        vm.set_variable("b", 3.0);
        vm.set_variable("c", 2.0);
        vm.set_variable("d", 4.0);
        vm.execute(prog);
        VarSlot res = vm.get_result();
        double val = std::get<double>(res.value);
        REQUIRE( std::abs(val - 1.0) < 1e-9 ); // true
    }
}


