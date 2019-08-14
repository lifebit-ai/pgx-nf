#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/pgx-nf
========================================================================================
 lifebit-ai/pgx-nf Pharamacogenomics Pipeline.
 #### Homepage / Documentation
 https://github.com/lifebit-ai/pgx-nf
----------------------------------------------------------------------------------------
*/

Channel
  .fromPath(params.tissue_annotation)
  .ifEmpty { exit 1, "Tissue annotation CSV file not found: ${params.tissue_annotation}" }
  .set { tissue_annotation }
// Channel
//   .fromPath(params.cellcuration)
//   .ifEmpty { exit 1, "Cell curation RData file not found: ${params.cellcuration}" }
//   .into { cellcuration }

/*--------------------------------------------------
  Compile cell curation
---------------------------------------------------*/

process compilecellcuration {

  tag "$tissue_annotation"
  container 'bhklab/pharmacogxcwl'

  input:
  file(tissue_annotation) from tissue_annotation
  // file(cellcuration) from cellcuration

  output:
  file("tissue_cur.RData") into tissuecuration

  script:
  """
  tissue_curation.R
  """
}