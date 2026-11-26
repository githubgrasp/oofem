#
# Restart test for interface bond laws and CDPM2: the final step after
# restarting from a context file must match the uninterrupted analysis.
#   bondceb02        - intelpoint + bondceb (interface status jump/traction)
#   linkslip02       - intelpoint + linkslip
#   con2dpm_random01 - CDPM2 with randvars (random properties in context)
#
OOFEM=$1
shift
OOFEM_ARGS=("$@")
cd "$(dirname "$0")" || exit 1
. ../restartcheck.sh

restart_case bondceb02.in 10
restart_case linkslip02.in 10
restart_case con2dpm_random01.in 4

exit $RESTART_FAILED
