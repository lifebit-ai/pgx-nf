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
Channel
  .fromPath(params.drug_annotation)
  .ifEmpty { exit 1, "Drug annotation CSV file not found: ${params.drug_annotation}" }
  .set { drug_annotation }
Channel
  .fromPath(params.published_info)
  .ifEmpty { exit 1, "Published info XLSX file not found: ${params.published_info}" }
  .set { published_info }
Channel
  .fromPath(params.drug_raw)
  .ifEmpty { exit 1, "Raw drug dose reponse CSV file not found: ${params.drug_raw}" }
  .set { drug_raw }
Channel
  .fromPath(params.drug_conc)
  .ifEmpty { exit 1, "Raw drug dose reponse CSV file not found: ${params.drug_conc}" }
  .set { drug_conc }
Channel
  .fromPath(params.gr_values)
  .ifEmpty { exit 1, "GR values TSV file not found: ${params.gr_values}" }
  .set { gr_values }
Channel
  .fromPath(params.crosscell)
  .ifEmpty { exit 1, "Cross referencing cells TXT file not found: ${params.crosscell}" }
  .set { crosscell }
Channel
  .fromPath(params.crossdrug)
  .ifEmpty { exit 1, "Cross referencing drug TXT file not found: ${params.crossdrug}" }
  .set { crossdrug }
Channel
  .fromPath(params.matrix)
  .ifEmpty { exit 1, "RNA-Seq matrix TXT file not found: ${params.matrix}" }
  .set { matrix }
Channel
  .fromPath(params.rnaseqinfo)
  .ifEmpty { exit 1, "RNA-Seq info CSV file not found: ${params.rnaseqinfo}" }
  .set { rnaseqinfo }
Channel
  .fromPath(params.rnaseqfeature)
  .ifEmpty { exit 1, "RNA-Seq feature CSV file not found: ${params.rnaseqfeature}" }
  .set { rnaseqfeature }
Channel
  .fromPath(params.expression)
  .ifEmpty { exit 1, "RNA-Seq expression TXT file not found: ${params.expression}" }
  .set { expression }
Channel
  .fromPath(params.counts)
  .ifEmpty { exit 1, "RNA-Seq counts TXT file not found: ${params.counts}" }
  .set { counts }
  

/*--------------------------------------------------
  Compile cell curation
---------------------------------------------------*/

process compilecellcuration {

  tag "$cell_annotation"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cell_annotation) from cell_annotation

  output:
  set file("cell_cur.RData") into cellcuration_tissue, cellcuration_cellline, cellcuration_recomput, cellcuration_rnaseq

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
  file(cellcuration) from cellcuration_tissue
  file(tissue_annotation) from tissue_annotation

  output:
  file("tissue_cur.RData") into tissuecuration

  script:
  """
  tissue_curation.R $tissue_annotation $cellcuration
  """
}

/*--------------------------------------------------
  Compile drug curation
---------------------------------------------------*/

process compiledrugcuration {

  tag "$drug_annotation"
  container 'bhklab/pharmacogxcwl'

  input:
  file(drug_annotation) from drug_annotation

  output:
  file("drug_cur.RData") into drugcuration

  script:
  """
  drug_curation.R $drug_annotation
  """
}

/*--------------------------------------------------
  Compile cell line info
---------------------------------------------------*/

process compilecelllineinfo {

  tag "$cellcuration"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cellcuration) from cellcuration_cellline
  file(published_info) from published_info

  output:
  file("cellline_info.RData") into celllineinfo

  script:
  """
  cellline_info.R $published_info $cellcuration
  """
}

/*--------------------------------------------------
  Recomputation
---------------------------------------------------*/

process recomputation {

  tag "${cellcuration},${drugcuration}"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cellcuration) from cellcuration_recomput
  file(drugcuration) from drugcuration
  file(drug_raw) from drug_raw
  file(drug_conc) from drug_conc
  file(gr_values) from gr_values
  file(crosscell) from crosscell
  file(crossdrug) from crossdrug

  output:
  file("drug_norm_post2017.RData") into drugnormpost
  file("GRAYrecomputed_2017.RData") into gray_recomputed2017

  script:
  """
  recomputed_2017.R \
    $cellcuration \
    $drugcuration \
    $drug_raw \
    $drug_conc \
    $gr_values \
    $crosscell \
    $crossdrug
  """
}

/*--------------------------------------------------
  Compile RNA-Seq
---------------------------------------------------*/

process compileRNAseq {

  tag "$cellcuration"
  container 'bhklab/pharmacogxcwl'

  input:
  file(cellcuration) from cellcuration_rnaseq
  file(matrix) from matrix
  file(rnaseqinfo) from rnaseqinfo
  file(rnaseqfeature) from rnaseqfeature
  file(expression) from expression
  file(counts) from counts

  output:
  file("RNAseq_processed.RData") into rnaseq

  script:
  """
  rna_seq.R \
    $cellcuration \
    $matrix \
    $rnaseqinfo \
    $rnaseqfeature \
    $expression \
    $counts
  """
}