#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include nf-core modules
include { TIARA_TIARA                                                  } from '../../modules/nf-core/tiara/tiara/main'
include { DASTOOL_FASTATOCONTIG2BIN as DASTOOL_FASTATOCONTIG2BIN_TIARA } from '../../modules/nf-core/dastool/fastatocontig2bin/main'
include { CSVTK_CONCAT as TIARA_SUMMARY                                } from '../../modules/nf-core/csvtk/concat/main'

// Include local modules
include { TIARA_CLASSIFY                                               } from '../../modules/local/tiara_classify/main'

workflow DOMAIN_CLASSIFICATION {
    take:
    ch_assemblies // channel: [ val(meta), path(assembly) ]
    ch_bins       // channel: [ val(meta), [path(bins)] ] - GROUPED bins as list

    main:
    ch_versions = Channel.empty()

    // Run Tiara on assemblies for domain classification
    TIARA_TIARA(ch_assemblies)
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions)

    // Add bin type metadata to distinguish bins
    ch_bins_with_type = ch_bins.map { meta, bin_list ->
        def meta_new = meta.clone()
        // Determine if these are refined or unrefined bins based on metadata
        meta_new.bin = meta.refinement == 'refined' ? 'refined_bins' : 'bins'
        // Ensure bin_list is a list (flatten if nested)
        def bins = bin_list instanceof List ? bin_list.flatten() : [bin_list]
        [meta_new, bins]
    }

    // Generate contig2bin files for each bin group
    DASTOOL_FASTATOCONTIG2BIN_TIARA(ch_bins_with_type, 'fa')
    ch_versions = ch_versions.mix(DASTOOL_FASTATOCONTIG2BIN_TIARA.out.versions)

    // FIX: Follow MAG pipeline pattern - combine bins with contig2bin FIRST
    ch_contigs_to_bin_tiara = DASTOOL_FASTATOCONTIG2BIN_TIARA.out.fastatocontig2bin
        .combine(ch_bins_with_type, by: 0)
        .map { meta, contig2bin, bin_list ->
            // Create meta_join by removing binner-specific fields
            def meta_join = meta.clone()
            meta_join.remove('binner')
            meta_join.remove('bin')
            meta_join.remove('refinement')
            meta_join.remove('bin_id')
            [meta_join, meta, contig2bin, bin_list]
        }

    // FIX: Combine with Tiara classifications
    ch_tiara_classify_input = ch_contigs_to_bin_tiara
        .combine(TIARA_TIARA.out.classifications, by: 0)
        .map { _meta_join, meta, contig2bin, bin_list, classifications ->
            // Return proper tuple for TIARA_CLASSIFY
            [meta, classifications, contig2bin, bin_list]
        }

    // Classify bins using Tiara results
    TIARA_CLASSIFY(ch_tiara_classify_input)
    ch_versions = ch_versions.mix(TIARA_CLASSIFY.out.versions)

    // Process classified outputs - need to transpose to get individual bins back
    ch_eukarya_bins = TIARA_CLASSIFY.out.eukarya_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'eukarya'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    ch_prokarya_bins = TIARA_CLASSIFY.out.prokarya_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'prokarya'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    ch_bacteria_bins = TIARA_CLASSIFY.out.bacteria_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'bacteria'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    ch_archaea_bins = TIARA_CLASSIFY.out.archaea_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'archaea'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    ch_organelle_bins = TIARA_CLASSIFY.out.organelle_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'organelle'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    ch_unknown_bins = TIARA_CLASSIFY.out.unknown_bins
        .transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.domain = 'unknown'
            meta_new.remove('bin')
            meta_new.bin_id = bin.baseName
            [meta_new, bin]
        }

    // Combine all classified bins
    ch_classified_bins = ch_eukarya_bins
        .mix(ch_prokarya_bins)
        .mix(ch_bacteria_bins)
        .mix(ch_archaea_bins)
        .mix(ch_organelle_bins)
        .mix(ch_unknown_bins)

    // FIX: Collect bin classifications properly for TIARA_SUMMARY
    // Extract just the classification files and create a proper tuple
    ch_bin_classifications = TIARA_CLASSIFY.out.bin_classifications
        .map { meta, classification -> classification }  // Extract just files
        .collect()  // Collect all files into a list
        .map { files -> [[:], files] }  // Create tuple with empty meta and file list

    TIARA_SUMMARY(ch_bin_classifications, 'tsv', 'tsv')
    ch_versions = ch_versions.mix(TIARA_SUMMARY.out.versions)

    emit:
    classified_bins = ch_classified_bins
    versions        = ch_versions
}
