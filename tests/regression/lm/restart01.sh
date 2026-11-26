#
# Restart test for lattice materials: the final step after restarting
# from a context file must match the uninterrupted analysis.
#   latticeplastdamvisco_mps_1 - LatticePlasticityDamage + viscoelastic slave
#   latticedamagevisco_mps_1   - LatticeDamage + viscoelastic slave
#   latticeviscoelastic_mps_1  - LatticeViscoelastic
#   latticeframesteelplastic1  - LatticeFrameSteelPlastic (hardening kappa)
#   latticelink3dnl_assembly   - latticelink3dnl + LatticeSlip
#
OOFEM=$1
shift
OOFEM_ARGS=("$@")
cd "$(dirname "$0")" || exit 1
. ../restartcheck.sh

restart_case latticeplastdamvisco_mps_1.in 12
restart_case latticedamagevisco_mps_1.in 12
restart_case latticeviscoelastic_mps_1.in 12
restart_case latticeframesteelplastic1.in 10
restart_case latticelink3dnl_assembly.in 8

exit $RESTART_FAILED
