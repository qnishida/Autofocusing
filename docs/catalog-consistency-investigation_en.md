# Catalog consistency and regression checks

Historical agreement and scientific correctness are separate checks. A catalog
may reproduce an earlier implementation that contains a known defect.

## Effect of the current correction

The spectral matrix is evaluated after fitting and does not feed detection.
Evaluating corrected and original self-term formulas at the same fitted points
therefore changes matrix columns 34–38 while retaining columns 1–33 by construction.
Independent end-to-end comparisons are still needed to detect changes elsewhere.

The original same-station subtraction was already present in the initial public
source. The adopted correction agrees with an explicit sum over distinct station
pairs; see [matrix verification](spectral-matrix-3c-2005-verification_en.md).
Agreement with an old catalog alone cannot validate the matrix formula.

## Restart history and QC

`eval_deri` retains process-global integration and gap state across windows,
segments and dates. Starting directly on a target day can produce different
accepted windows from a continuous run, even with identical waveform files.
A preceding-day control demonstrated this restart effect in a historical audit.
It does not establish that one warmup day always recovers a full archive's QC state.

Preserve the start date and preceding QC sequence when testing a regression.
Record accepted windows and initial candidates before comparing fitted outputs.

## Known earlier diagnostics fix

The original code left `prm_init` uninitialized. Commit `393a93e` initializes
it from the incoming fit. This affects fit-change diagnostics in columns 9–10,
not the optimization path. Uninitialized values must not be interpreted as
reproducible zeros or as evidence of a numerical-method change.

## Public-source control

A controlled comparison must identify source revisions, compatibility patches,
compiler/library versions, input hashes, FFTW plans and Bootstrap seeds.
Historical waveform hashes and complete compiler settings were unavailable in
the catalog investigation, so full input identity and a unique cause for every
remaining difference could not be established.

Use [catalog matching](historical-catalog-agreement_en.md) to report matched and
unmatched rows separately. Detailed event-level investigations and restart
experiments are maintained in the separate analysis workspace.
