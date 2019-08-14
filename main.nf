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
  .fromPath(params.cell_annotation)
  .ifEmpty { exit 1, "Cell annotation CSV file not found: ${params.cell_annotation}" }
  .set { cell_annotation }

/*--------------------------------------------------
  Compile cell curation
---------------------------------------------------*/

process compilecellcuration {

  tag "$cell_annotation"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cell_annotation) from cell_annotation

  output:
  file("cell_cur.RData") into cellcuration

  script:
  """
  cell_curation.R $cell_annotation
  """
}