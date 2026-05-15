# Schema compatibility matrix

| ngsPedigree Stage 3 version | Stage 3 schema | ngsTracts version | ngsTracts expects | compatible? |
| --------------------------- | -------------- | ----------------- | ----------------- | ----------- |
| 0.1                         | 0.1            | 0.1.x             | 0.1               | yes         |

## Compatibility rules

- **Patch version (0.1.0 → 0.1.x)**: backward compatible. ngsTracts of any
  0.1.x reads Stage 3 schema 0.1 outputs.
- **Minor version (0.1 → 0.2)**: additive only. New optional columns may
  appear in Stage 3 outputs; ngsTracts ignores them. Required columns
  do not change.
- **Major version (0.x → 1.0, 1.0 → 2.0)**: breaking. Column meaning or
  required set changes. ngsTracts refuses to run; you must upgrade or
  downgrade one side.

## How ngsTracts checks compatibility

On startup, `STEP_TRC_01_classify_intervals.py` reads `stage3.args.tsv`
and looks for the `schema_version` line. It:

1. Parses the version as `MAJOR.MINOR` (extra components ignored).
2. Compares against its own supported schema version set
   (currently `{"0.1"}`).
3. On mismatch: emits a clear error with both versions, the supported
   set, and the path to this doc. Refuses to run.
4. On match: proceeds. Logs the version to its own `.args` output.

## How to add a new column to the contract

If Stage 3 needs to emit something new that ngsTracts will read:

1. Bump Stage 3 schema_version (minor if additive, major if breaking).
2. Update `SCHEMA.md` to add the column.
3. Update this compatibility matrix.
4. If additive: ngsTracts can still consume the old schema without code
   change, but a later release can start recognizing the new column.
5. If breaking: bump ngsTracts major version, code the new reader,
   update the matrix.
