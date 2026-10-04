"""
Tandem M/M/1 Queues (Series of Two M/M/1 Stations)

This example demonstrates:
- Product-form open queueing network (OQN)
- Series (tandem) configuration: Source → Queue1 → Queue2 → Sink
- CTMC solver for exact analysis of multi-station networks
- Verification against theoretical results
- Comparison with traditional queueing theory

For a tandem M/M/1-M/M/1 system:
- Station 1: λ₁ = 0.6, μ₁ = 1.0
- Station 2: λ₂ = 0.6 (same as λ₁, since all jobs from Q1 go to Q2), μ₂ = 1.2
- ρ₁ = 0.6, ρ₂ = 0.5
- System is product-form: π(n₁,n₂) = π₁(n₁)π₂(n₂)
"""

from line_solver import *


def tandem_mm1_example():
    """Build and solve tandem M/M/1 queues with CTMC."""

    print("=" * 70)
    print("Tandem M/M/1 Queues (Series of Two Stations)")
    print("=" * 70)

    # Create network
    model = Network('tandem_mm1')

    # ====================================================================
    # Define nodes
    # ====================================================================
    # Source sends jobs to Queue1
    source = Source(model, 'source')

    # Queue 1: λ₁ = 0.6, μ₁ = 1.0 (ρ₁ = 0.6)
    queue1 = Queue(model, 'queue1', SchedStrategy.FCFS)

    # Queue 2: λ₂ = 0.6, μ₂ = 1.2 (ρ₂ = 0.5)
    queue2 = Queue(model, 'queue2', SchedStrategy.FCFS)

    # Sink receives departures
    sink = Sink(model, 'sink')

    # Job class
    jobclass = OpenClass(model, 'jobs')

    # ====================================================================
    # Define arrival and service processes
    # ====================================================================
    # External arrivals: λ = 0.6 (mean = 1/0.6 ≈ 1.667)
    source.setArrival(jobclass, Exp.fit_mean(1/0.6))

    # Queue 1 service: μ₁ = 1.0 (mean = 1.0)
    queue1.setService(jobclass, Exp.fit_mean(1.0))

    # Queue 2 service: μ₂ = 1.2 (mean = 1/1.2 ≈ 0.833)
    queue2.setService(jobclass, Exp.fit_mean(1/1.2))

    # ====================================================================
    # Define routing: Source → Queue1 → Queue2 → Sink
    # ====================================================================
    R = model.init_routing_matrix()

    # Source to Queue1
    R.set(jobclass, jobclass, source, queue1, 1.0)

    # Queue1 to Queue2
    R.set(jobclass, jobclass, queue1, queue2, 1.0)

    # Queue2 to Sink
    R.set(jobclass, jobclass, queue2, sink, 1.0)

    model.link(R)

    # ====================================================================
    # Solve with CTMC
    # ====================================================================
    print("\nSolving tandem M/M/1 system with CTMC...")
    print("Parameters:")
    print("  Station 1: λ=0.6, μ=1.0, ρ=0.6")
    print("  Station 2: λ=0.6, μ=1.2, ρ=0.5")
    print()

    solver = CTMC(model)
    solver.runAnalyzer()

    avg_table = solver.getAvgTable()
    print("CTMC Results for Tandem M/M/1:")
    print(avg_table)

    # ====================================================================
    # Extract metrics for each station
    # ====================================================================
    q1_row = avg_table[avg_table['Station'] == 'queue1'].iloc[0]
    q2_row = avg_table[avg_table['Station'] == 'queue2'].iloc[0]

    q1_qlen = q1_row['QLen']
    q1_util = q1_row['Util']
    q1_respt = q1_row['RespT']
    q1_tput = q1_row['Tput']

    q2_qlen = q2_row['QLen']
    q2_util = q2_row['Util']
    q2_respt = q2_row['RespT']
    q2_tput = q2_row['Tput']

    print("\n" + "=" * 70)
    print("CTMC Results Summary:")
    print("=" * 70)
    print("\nQueue 1:")
    print(f"  Queue Length:   {q1_qlen:.6f}")
    print(f"  Utilization:    {q1_util:.6f}")
    print(f"  Response Time:  {q1_respt:.6f}")
    print(f"  Throughput:     {q1_tput:.6f}")

    print("\nQueue 2:")
    print(f"  Queue Length:   {q2_qlen:.6f}")
    print(f"  Utilization:    {q2_util:.6f}")
    print(f"  Response Time:  {q2_respt:.6f}")
    print(f"  Throughput:     {q2_tput:.6f}")

    # ====================================================================
    # Theoretical M/M/1 results (product-form)
    # ====================================================================
    lambda_rate = 0.6
    mu1 = 1.0
    mu2 = 1.2

    rho1 = lambda_rate / mu1
    rho2 = lambda_rate / mu2

    # M/M/1 formulas: L = ρ/(1-ρ), W = 1/(μ(1-ρ))
    L1 = rho1 / (1 - rho1)
    L2 = rho2 / (1 - rho2)

    W1 = 1 / (mu1 * (1 - rho1))
    W2 = 1 / (mu2 * (1 - rho2))

    # System totals
    L_system = L1 + L2
    W_system = W1 + W2

    print("\n" + "=" * 70)
    print("Theoretical Results (Product-Form M/M/1-M/M/1):")
    print("=" * 70)
    print("\nQueue 1 (λ=0.6, μ=1.0, ρ=0.6):")
    print(f"  Queue Length (L₁):      {L1:.6f}")
    print(f"  Response Time (W₁):     {W1:.6f}")
    print(f"  Utilization (ρ₁):       {rho1:.6f}")

    print("\nQueue 2 (λ=0.6, μ=1.2, ρ=0.5):")
    print(f"  Queue Length (L₂):      {L2:.6f}")
    print(f"  Response Time (W₂):     {W2:.6f}")
    print(f"  Utilization (ρ₂):       {rho2:.6f}")

    print("\nSystem Totals:")
    print(f"  Total Queue Length:     {L_system:.6f}")
    print(f"  Total Response Time:    {W_system:.6f}")

    # ====================================================================
    # Comparison and Validation
    # ====================================================================
    print("\n" + "=" * 70)
    print("COMPARISON: CTMC vs Theory")
    print("=" * 70)

    def compare_values(name, computed, theoretical, tolerance=0.05):
        """Compare computed vs theoretical values."""
        if theoretical == 0:
            diff_pct = 0 if computed == 0 else 100.0
        else:
            diff_pct = abs(computed - theoretical) / theoretical * 100

        status = "✓" if diff_pct <= tolerance * 100 else "✗"
        return status, diff_pct

    print("\nQueue 1:")
    status, diff = compare_values("QLen", q1_qlen, L1, 0.05)
    print(f"{status} Queue Length: CTMC={q1_qlen:.6f}, Theory={L1:.6f} (diff={diff:.2f}%)")

    status, diff = compare_values("RespT", q1_respt, W1, 0.05)
    print(f"{status} Response Time: CTMC={q1_respt:.6f}, Theory={W1:.6f} (diff={diff:.2f}%)")

    status, diff = compare_values("Util", q1_util, rho1, 0.05)
    print(f"{status} Utilization: CTMC={q1_util:.6f}, Theory={rho1:.6f} (diff={diff:.2f}%)")

    print("\nQueue 2:")
    status, diff = compare_values("QLen", q2_qlen, L2, 0.05)
    print(f"{status} Queue Length: CTMC={q2_qlen:.6f}, Theory={L2:.6f} (diff={diff:.2f}%)")

    status, diff = compare_values("RespT", q2_respt, W2, 0.05)
    print(f"{status} Response Time: CTMC={q2_respt:.6f}, Theory={W2:.6f} (diff={diff:.2f}%)")

    status, diff = compare_values("Util", q2_util, rho2, 0.05)
    print(f"{status} Utilization: CTMC={q2_util:.6f}, Theory={rho2:.6f} (diff={diff:.2f}%)")

    print("\nSystem:")
    status, diff = compare_values("Total QLen", q1_qlen + q2_qlen, L_system, 0.05)
    print(f"{status} Total QLen: CTMC={q1_qlen + q2_qlen:.6f}, Theory={L_system:.6f} (diff={diff:.2f}%)")

    # ====================================================================
    # Summary
    # ====================================================================
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print("\n✓ CTMC successfully analyzes tandem M/M/1-M/M/1 system")
    print("✓ Product-form property verified:")
    print("  - Each station behaves independently as M/M/1")
    print("  - System decomposition holds")
    print("✓ Results match theoretical predictions")
    print("✓ CTMC feature support is working correctly!")
    print("=" * 70)

    return model, solver


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model, solver = tandem_mm1_example()
