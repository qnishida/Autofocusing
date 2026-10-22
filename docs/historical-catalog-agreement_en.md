# Historical catalog comparison

## Counts and matching

Match events by equal start time, end time and detected component. Within each
group, first maximize the number of one-to-one matches within a specified
Euclidean distance in `(px, py)`, then minimize the total slowness distance.
The supplied comparison tool defaults to 0.005 s/km, one search-grid step.
Power is not used in matching. This operational correspondence does not prove
physical event identity.

Report matched counts, unmatched rows and tolerance sensitivity. For matched
events, report slowness distance and relative differences with an explicit
denominator. Small medians do not establish agreement for every event.
Do not silently apply a slowness sign reversal or waveform-orbit equivalence.

## Input identity and other columns

Compare waveform and catalog hashes, run start/QC history, accepted windows,
component mode, frequency band, search grid, FFTW planning and Bootstrap seeds.
If historical hashes or settings are missing, input identity remains unverified.
A shared 38-column layout does not imply unchanged diagnostics or acceptance.

In particular, the `prm_init` fix changed fit-change diagnostics, and the
spectral-matrix correction changed post-fit matrix output. See the
[consistency checks](catalog-consistency-investigation_en.md).

## Reproduction

```bash
python3 tests/compare_historical_catalog.py --help
```

Supply the historical catalog, an expanded audit JSON and a new report path.
Both catalogs are external inputs; they are not distributed with the software.
Keep pair identifiers, unmatched rows and detailed output outside the public
source checkout. Historical event-level results are kept in the analysis workspace.
