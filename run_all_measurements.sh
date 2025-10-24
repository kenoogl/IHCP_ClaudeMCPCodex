#!/bin/bash

# Step 1 remaining tests (3 patterns)
echo "=== Step 1 Remaining (3 tests) ==="
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl --solver pbicgstab --precond gs --nt 10 --basesize 5000 2>&1 | tee shared/results/threads_basesize/step1_t4_bs5000.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/test_dhcp_solver.jl --solver pbicgstab --precond gs --nt 10 --basesize 2000 2>&1 | tee shared/results/threads_basesize/step1_t8_bs2000.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/test_dhcp_solver.jl --solver pbicgstab --precond gs --nt 10 --basesize 5000 2>&1 | tee shared/results/threads_basesize/step1_t8_bs5000.log

# Step 2 all tests (18 patterns: 3 threads × 6 basesizes)
echo "=== Step 2: CGM Full (18 tests) ==="

# 2 threads
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 400 2>&1 | tee shared/results/threads_basesize/step2_t2_bs400.log
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 500 2>&1 | tee shared/results/threads_basesize/step2_t2_bs500.log
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 600 2>&1 | tee shared/results/threads_basesize/step2_t2_bs600.log
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 1000 2>&1 | tee shared/results/threads_basesize/step2_t2_bs1000.log
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 2000 2>&1 | tee shared/results/threads_basesize/step2_t2_bs2000.log
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 5000 2>&1 | tee shared/results/threads_basesize/step2_t2_bs5000.log

# 4 threads
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 400 2>&1 | tee shared/results/threads_basesize/step2_t4_bs400.log
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 500 2>&1 | tee shared/results/threads_basesize/step2_t4_bs500.log
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 600 2>&1 | tee shared/results/threads_basesize/step2_t4_bs600.log
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 1000 2>&1 | tee shared/results/threads_basesize/step2_t4_bs1000.log
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 2000 2>&1 | tee shared/results/threads_basesize/step2_t4_bs2000.log
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 5000 2>&1 | tee shared/results/threads_basesize/step2_t4_bs5000.log

# 8 threads
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 400 2>&1 | tee shared/results/threads_basesize/step2_t8_bs400.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 500 2>&1 | tee shared/results/threads_basesize/step2_t8_bs500.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 600 2>&1 | tee shared/results/threads_basesize/step2_t8_bs600.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 1000 2>&1 | tee shared/results/threads_basesize/step2_t8_bs1000.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 2000 2>&1 | tee shared/results/threads_basesize/step2_t8_bs2000.log
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs --basesize 5000 2>&1 | tee shared/results/threads_basesize/step2_t8_bs5000.log

echo "=== All measurements completed! ==="
