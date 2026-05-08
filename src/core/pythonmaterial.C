#include "pythonmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gausspoint.h"
#include "timestep.h"
#include "error.h"
#include "classfactory.h"

#ifdef _USE_NANOBIND
    #include <nanobind/nanobind.h>
    #include "../../bindings/python/oofemarray-nanobind.h"
    namespace nb = nanobind;
#elif defined(_PYBIND_BINDINGS)
    #include <pybind11/embed.h>
    #include <pybind11/numpy.h>
    namespace py = pybind11;
#endif

namespace oofem {

REGISTER_Material(PythonMaterial);

PythonMaterialStatus::PythonMaterialStatus(GaussPoint * gp) : MaterialStatus(gp)
{
}

void PythonMaterialStatus::initTempStatus()
{
    MaterialStatus::initTempStatus();
#ifdef _USE_NANOBIND
    if (stateDict.ptr() != nullptr) {
        tempStateDict = nb::cast<nb::dict>(stateDict.attr("copy")());
    } else {
        tempStateDict = nb::dict();
    }
#elif defined(_PYBIND_BINDINGS)
    if (stateDict.ptr() != nullptr) {
        tempStateDict = py::cast<py::dict>(stateDict.attr("copy")());
    } else {
        tempStateDict = py::dict();
    }
#endif
}

void PythonMaterialStatus::updateYourself(TimeStep *tStep)
{
    MaterialStatus::updateYourself(tStep);
#ifdef _USE_NANOBIND
    if (tempStateDict.ptr() != nullptr) {
        stateDict = nb::cast<nb::dict>(tempStateDict.attr("copy")());
    } else {
        stateDict = nb::dict();
    }
#elif defined(_PYBIND_BINDINGS)
    if (tempStateDict.ptr() != nullptr) {
        stateDict = py::cast<py::dict>(tempStateDict.attr("copy")());
    } else {
        stateDict = py::dict();
    }
#endif
}

PythonMaterial::PythonMaterial(int n, Domain *d) : Material(n, d)
{
}

void PythonMaterial::initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    Material::initializeFrom(ir);
    IR_GIVE_FIELD(ir, moduleName, _IFT_PythonMaterial_module);
    IR_GIVE_FIELD(ir, objectName, _IFT_PythonMaterial_object);
}

std::unique_ptr<MaterialStatus> PythonMaterial::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<PythonMaterialStatus>(gp);
}

int PythonMaterial::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
    std::string key = std::to_string(type);

#ifdef _USE_NANOBIND
    if (ms->giveStateDictionary().contains(nb::str(key.c_str()))) {
        nb::object val = ms->giveStateDictionary()[nb::str(key.c_str())];
        try {
            answer.resize(1);
            answer.at(1) = nb::cast<double>(val);
            return 1;
        } catch (...) {}
        try {
            answer = nb::cast<FloatArray>(val);
            return 1;
        } catch (...) {}
        OOFEM_WARNING("Dictionary entry of material state not double or FloatArray");
    }
#elif defined(_PYBIND_BINDINGS)
    if (ms->giveStateDictionary().contains(key)) {
        py::object val = ms->giveStateDictionary()[key.c_str()];
        try {
            answer.resize(1);
            answer.at(1) = val.cast<double>();
            return 1;
        } catch (...) {}
        try {
            answer = val.cast<FloatArray>();
            return 1;
        } catch (...) {}
        OOFEM_WARNING("Dictionary entry of material state not double or FloatArray");
    }
#endif

#endif

    return Material::giveIPValue(answer, gp, type, tStep);
}

bool PythonMaterial::hasMaterialModeCapability(MaterialMode mode) const
{
#ifdef _USE_NANOBIND
    nb::module_ calc = nb::module_::import_(moduleName.c_str());
    nb::object result = calc.attr(objectName.c_str()).attr("hasMaterialModeCapability")(nb::cast(mode));
    return nb::cast<bool>(result);
#elif defined(_PYBIND_BINDINGS)
    py::module calc = py::module::import(moduleName.c_str());
    py::object result = calc.attr(objectName.c_str()).attr("hasMaterialModeCapability")(mode);
    return result.cast<bool>();
#else
    OOFEM_ERROR("Not compiled with python support.");
    return false;
#endif
}

void PythonMaterial::giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::module_ calc = nb::module_::import_(moduleName.c_str());
    nb::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicMatrix")(nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = nb::cast<FloatMatrix>(result);
#elif defined(_PYBIND_BINDINGS)
    py::module calc = py::module::import(moduleName.c_str());
    py::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicMatrix")(type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = result.cast<FloatMatrix>();
#else
    OOFEM_ERROR("Not compiled with python support.");
#endif
}

void PythonMaterial::giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::module_ calc = nb::module_::import_(moduleName.c_str());
    nb::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicVector")(nb::cast(flux), nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = nb::cast<FloatArray>(result);
#elif defined(_PYBIND_BINDINGS)
    py::module calc = py::module::import(moduleName.c_str());
    py::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicVector")(flux, type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    answer = result.cast<FloatArray>();
#else
    OOFEM_ERROR("Not compiled with python support.");
#endif
}

double PythonMaterial::giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const
{
#if defined(_USE_NANOBIND) || defined(_PYBIND_BINDINGS)
    auto ms = static_cast<PythonMaterialStatus *>(this->giveStatus(gp));
#endif

#ifdef _USE_NANOBIND
    nb::module_ calc = nb::module_::import_(moduleName.c_str());
    nb::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicValue")(nb::cast(type), nb::cast(gp), nb::cast(tStep), ms->giveStateDictionary(), ms->giveTempStateDictionary());
    return nb::cast<double>(result);
#elif defined(_PYBIND_BINDINGS)
    py::module calc = py::module::import(moduleName.c_str());
    py::object result = calc.attr(objectName.c_str()).attr("giveCharacteristicValue")(type, gp, tStep, ms->giveStateDictionary(), ms->giveTempStateDictionary());
    return result.cast<double>();
#else
    OOFEM_ERROR("Not compiled with python support.");
    return 0.0;
#endif
}

} // end namespace oofem
