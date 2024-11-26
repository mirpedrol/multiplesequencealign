include { UNTAR } from '../../../modules/nf-core/untar/main'

workflow EXTRACT_STRUCTURES {
    take:
    ch_structures // channel: [ meta, /path/to/file.pdb ]

    main:

    // Structures are taken from a directory of PDB files.
    // If the directory is compressed, it is uncompressed first.
    ch_structures
        .branch { structures ->
            compressed:   structures[1].name.endsWith('.tar.gz')
            uncompressed: true
        }
        .set { ch_structures_branched }

    UNTAR(ch_structures_branched.compressed)
        .untar
        .mix(ch_structures_branched.uncompressed)
        .map {
            meta,dir ->
                [ meta,file(dir).listFiles().collect() ]
        }
        .set { ch_structures_list }

    emit:
    structures = ch_structures_list
}
