# Compatibility policy

Versions follow Semantic Versioning. What a version number promises:

| surface | changes allowed in a patch | in a minor | only in a major |
| --- | --- | --- | --- |
| `cells.tsv` / `cells.json` column names and meaning | none | new columns appended | renaming, removal, change of definition |
| `spliceqc.tsv` pipeline contract | none | new columns appended, `contract_version` bumped | column removal or reorder |
| `summary.json`, `pipeline_step.json` keys | none | new keys | removal or type change |
| flag semantics (thresholds, FDR, direction) | none | none | any change |
| geneset catalog membership (`resources/genesets`) | none | genes or genesets added (catalog minor) | genes or genesets removed or renamed |
| stage-15 panel membership (`panel_version`) | none | none | any change, with a new `panel_version` |
| reference file format (`ref.json`) | none | new optional fields; older files stay readable | required fields, `version` bump |
| expression cache (`expr.bin`) | none | new format version with migration or rebuild | - |
| CLI flags | none | new flags, new defaults only when outputs are unchanged | removal, renamed flags, changed defaults |
| determinism (byte-identical outputs for the same inputs, flags and version) | guaranteed | guaranteed | guaranteed |

Before 1.0 a minor release may contain breaking changes; each one is listed
under "Changed (breaking)" in the changelog with the old and the new name.
The v0.3 rename of expression signatures to `*_expr` is the precedent.

A deprecation is announced one minor release before removal: the old name
keeps working and a warning names the replacement.
