# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Unit sheet

Add samples to `config/units.tsv`.
* Each unit has a `unit_name`. This can be a running number, or an actual run, lane or replicate id (for now this is not functional as we want to make this workflow to comply with the dna-seq-varlociraptor-workflow because we will import this wokrflow as a module there. For that reason, it's not called `samples.tsv`)
* Each unit has a `sample_name`, which associates it with the biological sample it comes from.
* For each unit, you need to specify:
  * `fq1` and `fq2` for paired end reads. These can point to any FASTQ files on your system.
 
