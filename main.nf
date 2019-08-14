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
  .into { cell_annotation; tissue_annotation }

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

/*--------------------------------------------------
  Compile tissue curation
---------------------------------------------------*/

process compiletissuecuration {

  tag "$cellcuration"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cellcuration) from cellcuration
  file(tissue_annotation) from tissue_annotation

  output:
  file("tissue_cur.RData") into tissuecuration

  script:
  """
  tissue_curation.R $tissue_annotation $cellcuration
  """
}
