# Spectral-matrix comparison after event selection

## Selection and interpretation

To isolate a matrix-formula change, evaluate both formulas using identical
waveforms, fitted parameters, weights, selected events and QC history.
Apply event selection to the detected component's score, with the threshold
and strictness stated explicitly. A historical audit used strict max/MAD > 7
for R/T and > 35 for U; these are audit choices, not new production defaults.

For a power ratio, include only events with positive numerator and denominator
in both formulas. State the ratio direction, denominator, excluded counts and
component groups. Corrected self-term-subtracted matrices need not be positive
semidefinite for finite data; negative powers must not be silently clipped.

## Test method

`tests/spectral_matrix.cpp` checks the production formula against explicit
ordered pairs of distinct stations, including Hermitian structure and the
diagonal/scalar relationship. See the
[matrix verification summary](spectral-matrix-3c-2005-verification_en.md).

`tests/run_matrix_audit_events.py` and `tests/analyze_matrix_ratios.py` support
real-data audits with separately obtained waveform and catalog inputs. Use
fixed test seeds and FFTW planning, and verify non-matrix output columns,
accepted windows and initial candidates before attributing differences to the matrix.

Historical-catalog agreement additionally requires the
[restart and input-identity controls](catalog-consistency-investigation_en.md).
Use consecutive-date controls as well as catalog-guided dates; do not count
overlapping reruns as independent samples.

## Reproduction

```bash
python3 tests/run_matrix_audit_events.py --help
python3 tests/analyze_matrix_ratios.py --help
```

Save event-level audits, input hashes and comparisons in an external result
directory. Detailed date selections and catalog statistics are maintained in
the separate analysis workspace. For statistical interpretation, see
[estimate distributions](matrix-estimate-distributions_en.md).
