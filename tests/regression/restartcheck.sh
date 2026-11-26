#
# Helpers for restart regression tests (sourced by */restart*.sh).
#
# restart_case <input> <restart step>
#   1. runs <input> uninterrupted with -c (context written every step);
#      the errorcheck block of the input is applied as usual,
#   2. restarts from <restart step> with -r; errorcheck rules for steps
#      that are no longer computed are removed,
#   3. compares the final step block of both .out files numerically.
# Requires OOFEM (executable) and OOFEM_ARGS (array) to be set.
#

RESTART_RTOL=${RESTART_RTOL:-1.e-3}
RESTART_FAILED=0

# prints the output block of the last step of an .out file
restart_last_block() {
    awk '
        /^Output for time/ { n = 0; inblock = 1 }
        /^User time consumed by solution step/ { inblock = 0 }
        # iteration counts may legitimately differ, results may not
        /iterations/ { next }
        inblock { buf[++n] = $0 }
        END { for ( i = 1; i <= n; i++ ) print buf[i] }
    ' "$1"
}

# compares two blocks token by token, numbers with relative tolerance
restart_compare() {
    awk -v rtol="$RESTART_RTOL" '
        function abs(v) { return v < 0 ? -v : v }
        function isnum(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/ }
        NR == FNR { ref[FNR] = $0; nref = FNR; next }
        {
            nrst = FNR
            if ( FNR > nref ) { print "  extra line in restart: " $0; bad = 1; exit }
            n1 = split(ref[FNR], x); n2 = split($0, y)
            ok = ( n1 == n2 )
            # noise floor relative to the largest value on the line
            scale = 0.
            for ( i = 1; i <= n1; i++ ) if ( isnum(x[i]) && abs(x[i]) > scale ) scale = abs(x[i])
            for ( i = 1; ok && i <= n1; i++ ) {
                if ( x[i] == y[i] ) continue
                if ( isnum(x[i]) && isnum(y[i]) && abs(x[i] - y[i]) <= rtol * abs(x[i]) + 1.e-6 * scale ) continue
                ok = 0
            }
            if ( !ok ) {
                if ( ++ndiff <= 10 ) print "  ref: " ref[FNR] "\n  rst: " $0
                bad = 1
            }
        }
        END {
            if ( nrst != nref ) { print "  block length differs: " nref " vs " nrst; bad = 1 }
            if ( nref == 0 ) { print "  empty reference block"; bad = 1 }
            if ( ndiff > 10 ) print "  ... " ndiff " differing lines in total"
            exit bad
        }
    ' "$1" "$2"
}

restart_case() {
    local input=$1 k=$2
    local name work out
    name=$(basename "$input" .in)
    work=$(mktemp -d "${TMPDIR:-/tmp}/oofem-restart-$name.XXXXXX")
    out=$(head -n 1 "$input" | tr -d '\r' | awk '{ print $1 }')

    echo "=== $input: restart from step $k (work dir $work)"

    # print every node and element so that the whole state is compared;
    # drop contextOutputStep, it overrides -c
    awk '
        tolower($1) == "outputmanager" { $0 = "OutputManager tstep_all dofman_all element_all" }
        {
            for ( i = 1; i < NF; i++ ) {
                if ( tolower($i) == "contextoutputstep" ) {
                    $i = ""; $(i + 1) = ""
                }
            }
            print
        }
    ' "$input" > "$work/full.in"

    # drop errorcheck rules of steps that the restarted run does not compute
    awk -v k="$k" '
        /^#%BEGIN_CHECK%/ { inchk = 1 }
        /^#%END_CHECK%/ { inchk = 0 }
        inchk { for ( i = 1; i < NF; i++ ) if ( $i == "tStep" && $(i + 1) + 0 <= k ) next }
        { print }
    ' "$work/full.in" > "$work/restart.in"

    if ! ( cd "$work" && "$OOFEM" "${OOFEM_ARGS[@]}" -f full.in -c > full.log 2>&1 ); then
        echo "  FAILED: uninterrupted run (log $work/full.log)"
        tail -n 20 "$work/full.log"
        RESTART_FAILED=1
        return
    fi
    restart_last_block "$work/$out" > "$work/ref.block"

    if ! ( cd "$work" && "$OOFEM" "${OOFEM_ARGS[@]}" -f restart.in -r "$k" > restart.log 2>&1 ); then
        echo "  FAILED: restarted run (log $work/restart.log)"
        tail -n 20 "$work/restart.log"
        RESTART_FAILED=1
        return
    fi
    restart_last_block "$work/$out" > "$work/rst.block"

    if restart_compare "$work/ref.block" "$work/rst.block"; then
        echo "  OK"
        rm -rf "$work"
    else
        echo "  FAILED: final step differs after restart (files kept in $work)"
        RESTART_FAILED=1
    fi
}
