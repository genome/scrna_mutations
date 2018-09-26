# Gene expression depth position tool

This tool takes a annotation file in genePred format, a Cell Ranger BAM file, and a list of cell barcodes and produces a file reporting on the base-pair position of alignments seen relative to the transcription start/stop sites. This tool was used to generate figure 1b in *Petti et al*.

## Output
After execution, the `outs` folder will contain a single CSV file with three columns -- `tx_id`, `tx_position`, and `read_molecule_fraction`. For each `tx_id`, there will be a `tx_position` for each transcript-coordinate position on the transcript, with the associated `read_molecule_fraction`. The `read_molecule_fraction` is calculated as the fraction of the number of unique molecules (cell-barcode/UMI pairs) seen at position X divided by the total seen across the entire transcript.

## Installation

In order to run this pipeline, you need to have [Martian](https://github.com/martian-lang/martian) and [CAT](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit) installed. Martian provides the pipeline execution environment, and CAT provides library modules that are used to parse the input transcript set.

## Execution

Modify the `example.mro` file provided to point to the outputs of your Cell Ranger run. The transcripts file must be in [genePred extended](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) format. You can use [UCSC tools](http://hgdownload.soe.ucsc.edu/admin/exe/) to convert to this format, including `gff3ToGenePred`. The genePred used in this study was extracted directly from the UCSC MySQL database with the command

```
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -Ne 'select * from wgEncodeGencodeBasicV27' hg38 | cut -f 2- > basic.gp
```

