# Changelog

## v2.0.0

- Added optional read-based post-processing modes (`all`, `phyloseq`, `taxburst`, and `rgi_bwt`) for both batch and single-sketch workflows, with prefixed outputs under `Readbased_Analysis/daa`.
- Added safe `postprocess_options` forwarding, explicit cross-parameter validation, and support for both `_1/_2` and `_R1/_R2` directory input pairing.
- Prepared configuration for Nextflow v26.04 strict syntax parsing.
- Replaced the custom `check_max` config helper with native Nextflow `resourceLimits`.
- Added default read-based publish directory parameters to avoid undefined-parameter warnings.
- Fixed duplicate input variable declarations in the read-based finalization module.
- Updated README, changelog, citation, and manifest metadata for the v2.0.0 / Nextflow v26.04 release target.
- Documented `conda`, `docker`, `apptainer`, and opt-in GPU profile usage without changing process containers.

Existing read-based and assembly output semantics remain unchanged when read-based post-processing is disabled.
