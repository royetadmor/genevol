#!/bin/bash
# Integration test suite for genevol.
# Usage (from project root): script/tests/run_tests.sh [--skip-sim]
#   --skip-sim  skip simulation phase and use existing data in sim_output/
set -e

PASS=0
FAIL=0
SIM_OUTPUT="./sim_output"
LOG_DIR="./log"
SKIP_SIM=false

for arg in "$@"; do
    [[ "$arg" == "--skip-sim" ]] && SKIP_SIM=true
done

mkdir -p "$LOG_DIR"

# Run one genevol test. Args: <display-name> <config-file> [<sim-output-dir>]
run_test() {
    local name="$1"
    local config="$2"
    local data_dir="${3:-$SIM_OUTPUT}"
    local log_file="$LOG_DIR/$(echo "$name" | tr ' /' '__').log"
    echo -n "  [$name] ... "
    docker run --rm \
        -v "$(realpath "$data_dir"):/app/genevol/sim_output:ro" \
        -v "$(realpath "$config"):/app/genevol/test_config.txt:ro" \
        genevol \
        ./GenEvol/genEvol "param=/app/genevol/test_config.txt" > "$log_file" 2>&1 || true

    if grep -q "GenEvol's done. Bye." "$log_file"; then
        echo "PASS"
        PASS=$((PASS + 1))
    else
        echo "FAIL  (see $log_file)"
        tail -30 "$log_file"
        FAIL=$((FAIL + 1))
    fi
}

echo "=== genevol integration tests ==="
echo ""

# ── Phase 0: Build ────────────────────────────────────────────────────────────
echo "[ Building ]"
echo -n "  [docker build] ... "
if docker build -t genevol -f Dockerfile . > "$LOG_DIR/build.log" 2>&1; then
    echo "PASS"
else
    echo "FAIL  (see $LOG_DIR/build.log)"
    tail -30 "$LOG_DIR/build.log"
    exit 1
fi
echo ""

# ── Phase 1: Generate test data ──────────────────────────────────────────────
echo "[ Generating test data ]"

if $SKIP_SIM; then
    echo "  (skipped — using existing data)"
else
    mkdir -p "$SIM_OUTPUT"

    echo -n "  [Simulate base] ... "
    if script/simulation/docker_run_simulator.sh test_configs/simulate_base.txt "$SIM_OUTPUT" \
            > "$LOG_DIR/simulate_base.log" 2>&1; then
        echo "PASS"
    else
        echo "FAIL  (see $LOG_DIR/simulate_base.log)"
        cat "$LOG_DIR/simulate_base.log"
        exit 1
    fi

    echo -n "  [Simulate WGD]  ... "
    if script/simulation/docker_run_simulator.sh test_configs/simulate_wgd.txt "$SIM_OUTPUT" \
            > "$LOG_DIR/simulate_wgd.log" 2>&1; then
        echo "PASS"
    else
        echo "FAIL  (see $LOG_DIR/simulate_wgd.log)"
        cat "$LOG_DIR/simulate_wgd.log"
        exit 1
    fi
fi

# ── Phase 2: Model fitting (against base simulation) ─────────────────────────
echo ""
echo "[ Model fitting ]"
run_test "CONST gain/loss"    "test_configs/const_model.txt"
run_test "LINEAR gain/loss"   "test_configs/linear_model.txt"
run_test "EXP gain/loss"      "test_configs/exp_model.txt"
run_test "Gamma rates"        "test_configs/gamma_rates.txt"
run_test "NegBinomial root"   "test_configs/negbinom_root.txt"
run_test "Large state space"  "test_configs/large_state_space.txt"

# ── Phase 3: WGD ─────────────────────────────────────────────────────────────
echo ""
echo "[ WGD ]"
run_test "WGD detect"         "test_configs/wgd_detect.txt"
run_test "WGD test (fixed q)" "test_configs/wgd_test_fixed.txt"
run_test "WGD test (free q)"  "test_configs/wgd_test_free.txt"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
if [ "$FAIL" -eq 0 ]; then
    echo "=== All $PASS tests passed ==="
else
    echo "=== Results: $PASS passed, $FAIL failed ==="
fi
[ "$FAIL" -eq 0 ]
