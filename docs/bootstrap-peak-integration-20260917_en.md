# Bootstrap/peak-search integration validation

The local peak-selection, double-precision and roundoff-termination changes were
committed as `dc38340`. Upstream added OpenMP three-component CPU
bootstrap in `49dbb27`. Running `git pull --no-rebase --no-edit origin main` on
`test/peak-search-boundaries` merged them without conflict as `b0ce192`.
The integration checks below compare the merged code with the final serial
peak/Newton version. They do not use the older, pre-fix catalog as the baseline.

The merged peak-search, Newton, roundoff and covariance-conversion routines are
text-identical to `dc38340`. The parallel bootstrap-batch routine is identical
to upstream. Parallelism evaluates independent samples, then reduces powers
in the original sample order. It does not change peak selection or 4×4 formulas.

The production executable and affected tests rebuild successfully. Thirteen
calculation/workspace tests pass, including the expanded bootstrap batch test.
The qualification-tools test also passes when rerun outside the sandbox:
its initial failure came from macOS `/usr/bin/time -l` calling a blocked
`sysctl kern.clockrate`, reproduced even with `/usr/bin/true`.
Thus 14 selected CTests have passing results. Raw evidence stays under the
ignored `docs/benchmarks/bootstrap-peak-integration-20260917/` directory.

## Real-data equivalence

The merged audit probe replayed Hi-net March 15–17, 2005, using three components,
0.1–0.25 Hz, slowness limits ±0.165 s/km, step 0.005 s/km, and eight CPU threads.
Waveform and catalog hashes, QC history, bootstrap seeds
(`1837 + 104729 * replicate`) and FFTW settings (`ESTIMATE | UNALIGNED`) matched
the saved final peak/Newton run. The full interval was replayed to preserve QC
history from the first day.

- All 145 event rows are byte-for-byte identical in all 38 columns, including
  component energies, covariance and bootstrap outputs; all output values are
  finite. Energy ratios derived from these outputs also match.
- All 147 initial candidates and all recorded optimizer steps, results and
  roundoff-convergence records are identical. Accepted windows also match.
- Selection is unchanged: 75 events (R31/T10/U34), using max/MAD > 7 for R/T and
  > 35 for U. The previously recovered R candidate is retained.
- Both the corrected spectral-matrix catalog and the audit's original-matrix
  control match their respective independent serial runs exactly.

These results support combining the parallel bootstrap with the validated peak
and precision changes. They do not establish validity for every possible event
or qualify a GPU backend.

## CPU timing

The timing comparison uses clean production sources at `dc38340` and `b0ce192`,
with only fixed bootstrap seeds and a shared deterministic FFTW object. The
extra numerical-audit calculations and convergence logging are excluded.
Both executables use the same Release/O3 Clang build and profiling settings.

On an Apple M4 Max (16 cores), March 15, 2005 was replayed with eight OpenMP
threads, three components and the same frequency/slowness settings. Each binary
received one excluded warmup followed by five measured runs in alternating
order, each in a fresh process. The OS cache was not purged, and no other agent
builds or tests ran concurrently. All twelve runs, including warmups, produced
the same 46 events byte for byte. Input hashes also remained unchanged.

| Median elapsed time | Before bootstrap parallelism | Combined version |
| --- | ---: | ---: |
| Whole process | 110.92 s | 56.14 s |
| Bootstrap stage | 64.73 s | 10.05 s |
| Grid stage | 4.149 s | 4.131 s |

Bootstrap was **6.44× faster**; the whole process took **49.39% less time**
(1.98× faster). Whole-process ranges were 110.54–117.23 s before and
55.88–56.17 s after; peak RSS was approximately 2.39 GiB for both.
All predefined gates passed: whole-process and grid medians within +5% of
baseline, at least 10% bootstrap improvement, and byte-identical event output.
These are measurements of this one-day CPU workload on this host. Summed
per-worker call times are not elapsed times and must not be used as speedups.

## Reproduction and scope

The validated merge is suitable for integration into main. Updating an installed
runtime executable requires the normal build/install procedure separately.

The local evidence directory contains the probe builder, comparison and timing
drivers, input/source/binary hashes, test logs and run reports. Detailed logs
and event catalogs remain Git-ignored. Historical probe builders depend on
specific source snapshots; rebuilding them directly against merged code is
not a substitute for the recorded integration driver.

Production bootstrap seeds still use `time(0) * i`. Faster execution can change
the sampled seeds even though the expression and resampling order are unchanged.
Use fixed seeds for correctness comparisons rather than interpreting ordinary
run-to-run bootstrap differences as a merge regression.
