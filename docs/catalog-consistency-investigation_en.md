# Catalog consistency and regression investigation — 2026-09-16

Detailed run logs, output catalogs and per-event comparison records are retained
locally and are not distributed in this repository. References to these records
below describe local artifacts.

[Japanese](catalog-consistency-investigation_jp.md)

The user identifies the historical catalog as a CPU result from essentially
the public-release source. Distinguish three objects: that historical catalog,
current CPU output before the self-term correction, and current CPU output
after the correction. Agreement with a historical implementation and correct
same-station subtraction are different validation questions.

## Effect of the current correction

Across all 677 current detections from 48 days, columns 1–33 are exactly equal
before and after; matrix columns 34–38 change in all 677 rows. Thus historical
matching, fitted parameters, max/MAD, window counts and bootstrap fields have
the same agreement with the historical catalog before and after this correction.
The audit evaluates both formulas on the same fit, so this equality is by
construction. The earlier independent original/corrected three-day runs also
matched all first 33 columns for their 45 events. Source inspection confirms
that matrix evaluation occurs after fitting and does not feed detection.

For fixed historical/current matched selected events, the median absolute
diagonal-power discrepancy is `100 * abs(current/historical - 1)`:

| Detected component / diagonal | N | Before correction | After correction |
| --- | ---: | ---: | ---: |
| R / SRR | 164 | 0.0056% | 4.10% |
| T / STT | 35 | 0.0079% | 10.04% |
| U / SUU | 149 | 0.0036% | 2.51% |

The corrected matrices depart further from the historical values. The incorrect
subtraction from `(icmp1, icmp1)` inside the `icmp2` loop is already present in
the initial public commit `e551a69`, not newly introduced by later optimization.
The proposed correction instead agrees with an independently calculated sum
over distinct station pairs. Historical agreement alone is not a correctness
criterion for these matrix entries.

Three-catalog numbers and row references
and [independent matrix verification](spectral-matrix-3c-2005-verification_en.md).

## Restart history explains part of the catalog discrepancy

`eval_deri` uses process-global `integ_pre`, `integ_old` and `count_gap`.
Initially `integ_pre=0`, making the first evaluated window fail the derivative
check. State persists across segments and dates. This behavior is present in
the initial public source. A run starting on a selected date does not have the
same QC history as a continuous archive run, even with identical waveforms.

All 46 matched historical/current rows with a one-window count difference
occur in the midnight segment of a current run's starting date (11 dates).
Thirty of those rows pass the detected-component selection. The historical
count is always one larger. The earlier catalog comparison did not account
for this restart condition; interpreting all differences as numerical rounding
or a code regression would be premature.

The same current audit executable was rerun March 15–17, allowing March 15 to
precede the target March 16–17. No algorithm or setting changed. On the 17
common historical matches on March 16:

| Check | Start March 16 | Start March 15 |
| --- | ---: | ---: |
| First-segment accepted windows | 47 | 48 (historical: 48) |
| Rows with historical window count | 13/17 | 17/17 |
| Median absolute max/MAD discrepancy | 0.1161% | 0.0926% |
| Maximum absolute max/MAD discrepancy | 1.1426% | 0.2135% |
| Median slowness-vector difference, s/km | 6.99e-7 | 4.21e-7 |

On March 17, the 20 common comparisons have unchanged statistics. This
demonstrates a restart effect, not full reproduction of the historical catalog
or a guarantee that one warm-up day always recovers a long run's QC state.
These overlapping reruns are excluded from the 48-day statistical sample.

Continuity comparison,
run checks and hashes.

## Known earlier diagnostics fix

Public code declared `prm_init` without initializing it. Commit `393a93e`
initializes it from the incoming fit, as documented at that time in
[horizontal verification](horizontal-verification.md). Only the fit-change
diagnostics use this variable after the optimization. Historical columns 9–10
behaving as though the initial p and distance were zero are consistent with
this defect; reading an uninitialized variable does not guarantee zero.
This is a known earlier bug fix, distinct from the matrix correction.

## Public-source control

The `src` trees at `149995f` (public release) and `e551a69` are identical.
A test binary was built from all public source files using the current LLVM
compiler and libraries. Compatibility changes are explicit includes, standard
namespace qualification and the macOS prefetch guard. FFTW plans and bootstrap
seeds are fixed as in the current audit. The known `prm_init` fix is backported
to avoid undefined diagnostic reads; no optimization formula is changed.

For March 16's first segment, four events match the current pre-correction
output byte for byte across all 38 columns. The public test driver stops after
that segment; normal segment length and calculations are retained.

Public control artifacts and compatibility patch.

The full November 4–5 public-source run also matches all 42 events byte for
byte across all 38 columns, using four threads as in the saved current run.
This includes the November 5 event whose slowness differs from the historical
catalog. Input waveform and source-catalog hashes match the current run and
remain unchanged. No public-to-current regression is reproduced in these
controlled examples; this does not qualify every day or execution mode.

Two-day public/current report,
public events,
current-before events,
complete test-driver patch.

Historical compiler flags, FFTW plans, random seeds and waveform hashes remain
unknown. Residual differences are not assigned to a specific cause without
evidence. The checks distinguish the matrix change from pre-existing catalog
differences; they do not establish absence of every possible regression.
