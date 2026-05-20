# Changelog

## v2.0.0

- Prepared configuration for Nextflow v26.04 strict syntax parsing.
- Replaced the custom `check_max` config helper with native Nextflow `resourceLimits`.
- Added default read-based publish directory parameters to avoid undefined-parameter warnings.
- Fixed duplicate input variable declarations in the read-based finalization module.
- Updated README, changelog, citation, and manifest metadata for the v2.0.0 / Nextflow v26.04 release target.
- Documented `conda`, `docker`, `apptainer`, and opt-in GPU profile usage without changing process containers.

No scientific logic, output semantics, or container references were intentionally changed.
