#ifndef pythonmaterial_h
#define pythonmaterial_h

#include "material.h"
#include "matstatus.h"
#include <string>
#include <memory>

#ifdef _USE_NANOBIND
    #include <nanobind/nanobind.h>
    namespace nb = nanobind;
#elif defined(_PYBIND_BINDINGS)
    #include <pybind11/pybind11.h>
    namespace py = pybind11;
#endif

#define _IFT_PythonMaterial_Name "pythonmaterial"
#define _IFT_PythonMaterial_module "module"
#define _IFT_PythonMaterial_object "object"

namespace oofem {

class OOFEM_EXPORT PythonMaterialStatus : public MaterialStatus
{
protected:
#ifdef _USE_NANOBIND
    nb::dict stateDict, tempStateDict;
#elif defined(_PYBIND_BINDINGS)
    py::dict stateDict, tempStateDict;
#endif

public:
    PythonMaterialStatus(GaussPoint * gp);
    ~PythonMaterialStatus() override = default;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

#ifdef _USE_NANOBIND
    nb::dict giveStateDictionary() { return stateDict; }
    nb::dict giveTempStateDictionary() { return tempStateDict; }
#elif defined(_PYBIND_BINDINGS)
    py::dict giveStateDictionary() { return stateDict; }
    py::dict giveTempStateDictionary() { return tempStateDict; }
#endif

    const char *giveClassName() const override { return "PythonMaterialStatus"; }
};

class OOFEM_EXPORT PythonMaterial : public Material
{
protected:
    std::string moduleName;
    std::string objectName;

public:
    PythonMaterial(int n, Domain *d);
    ~PythonMaterial() override = default;

    const char *giveClassName() const override { return "PythonMaterial"; }
    const char *giveInputRecordName() const override { return "PythonMaterial"; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    
    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;
    void giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;
    double giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const override;
};

} // end namespace oofem
#endif // pythonmaterial_h
