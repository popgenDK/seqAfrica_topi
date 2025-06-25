# a1kg genotype calling pipeline

An overview of how to use the genotype calling pipeline developed for the a1kg project.

1. [Install](#Install)
2. [Quickstart](#Quickstart)
3. [Input TSV](#Input-TSV)  
4. [Configuration](#Configuration)
5. [Targets](#Targets)
	1. [Filters](#Filters) 
	2. [Masks](#Masks) 
	3. [Merging](#Merging) 
	4. [Plink](#Plink) 
6. [SFS](#SFS)
	1. [Folding](#Folding) 
	2. [Plotting](#Plotting) 
	3. [Statistics](#Statistics) 
		1. [Fst](#Fst) 


## Quickstart
The pipeline uses standard Snakemake logic for determining the target from the CLI, and basic familiarity with Snakemake is recommended. The [docs](https://snakemake.readthedocs.io/en/stable/) are excellent. It may also be necessary to look at the source rules in `rules/*.smk` to find the name pattern of the targets to create.

If you're copying the code from elsewhere, ensure that you have the entire folder `gt/` and that this is the current working directory. The `snakefile` itself is *not* sufficient, since the logic is split across multiple files in `rules/` files, some rules depend on scripts in the `scripts/` subdirectory, and so on. 

Create an [input TSV](#Input-TSV) and a [configuration file](#Configuration), and you can run

```sh
snakemake path/to/target --configfile path/to/configfile -n
```

For a dry run, or 

```sh
nice snakemake path/to/target --configfile path/to/configfile -j $threads
```

to run with `$threads` threads using a default nice value. See [Targets](#Targets) for details on available targets.

## Input TSV

To point the pipeline to the input data, please construct a two-column TSV, in which the first column contains group names, and the second column contains paths to BAM files. For instance,

```
groupA	/path/to/individual1.bam
groupA	/path/to/individual2.bam
groupB	/path/to/individual3.bam
groupB	/path/to/individual4.bam
groupB	/path/to/individual5.bam
```

The group names serve two purposes. 

- First, all groups are called separately in the initial genotyping. That is, for the above, individuals 1 and 2 would be called together, and individuals 3, 4, and 5 would be called together.
- Second, one group can be defined as the "ingroup" in the `config.yaml` (see [Configuration](#Configuration)), which will define a special group for the purposes of later merging groups past genotype calling. See also [`Merging`](#Merging) below.

In this case, there is a species under study , i.e. topi, to be one group (and also the "ingroup"), and then various other groups for the outgroups (e.g. hartebeest, blue wildebeest, black wildebeest). It is perfectly possible to run with just a single group, if no split genotyping is preferred.

## Configuration

A YAML configuration file is required to run the pipeline. A template is provided under `configs/template.yaml`. Some notes on the keys:

- `ingroup`: Used for merging, must refer to a group name in the input TSV file. See [Input TSV](#Input-TSV) for more.
- `results_subdir`: Name of subdirectory under `results/` where generated files will appear.
- `ref`: Name of reference file (used for file name prefixes), as well as paths to reference and reference index.
- `sites` Path to sites file used for filtering, can be left blank for initial calling and added later for downstream sites filtering. The sites file is a BED file containing *good* sites, i.e. those that should be kept.
- `tsv`: Path to input TSV file. See [Input TSV](#Input-TSV).
- `populations`: A two-column TSV file linking sample names to populations. Only required for SFS, and can left blank. See (SFS)[#SFS] section for details.
- `tmp_dir`: Path to location for `tmp` files. Could be simply `tmp`.
  
*Important*: This run-specific configuration file is used in conjunction with default settings defined in the top-level `config.yaml`. Check that you agree with the cutoffs defined herein. To override these defaults, simply replicate the keys in your specific config file with the desired values. See examples in the existing YAML files.

## Targets
Some general pointers on targets to get started, however. Suppose we have set `results_subdir: animal` Then all files will be created in the `results/animal` sub-directory. In fact, with a few exceptions (plink files, in particular), most files will be created in `results/animal/vcf`. The main genotype call file is of the form `{refname}_{group}.bcf.gz`, so to call the `A` group (from the first column of the [input TSV](#Input-TSV)), the corresponding target will be `results/animal/vcf/Animal_A.bcf.gz` (where `Animal` is the mapping reference name, again from the input TSV). Then, to call all groups `A`, `B`, `C`, we could use a [brace expand](https://www.gnu.org/software/bash/manual/html_node/Brace-Expansion.html) and do

```sh
nice snakemake results/Goat/vcf/Goat_topi.bcf.gz --configfile path/to/configfile -j $threads
```

For convenience, common targets are defined the main snakefile using the [pseudo-rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#targets-and-aggregation) pattern. That is, we could have also run,

```sh
nice snakemake call_all_groups --configfile path/to/configfile -j $threads
```

to get the same result.


### Filters

Once you have a variant file, various filters and masks can be applied to your variants simply by adding suffixes to the file name. For instance, having generated `results/animal/vcf/Animal_A.bcf.gz` above, we could target `results/animal/vcf/Animal_A_variable.bcf.gz` to get a BCF file with only the variable sites. A (not necessarily exhaustive) list of current filter suffixes:

- `variable`: Remove monomorphic sites.
- `sites`: Remove sites not contained in sites BED file.
- `nomultiallelics`: Remove multiallelic sites.
- `noindels`: Remove indels.
- `missing`: Remove sites with any missing genotypes.

These can be composed in any order, e.g. `results/animal/vcf/Animal_A_variable_sites_missing.bcf.gz` will contain all variable sites after using the sites BED file with no missing genotypes.

Note that to maintain flexibility, the pipeline has to generate an intermediate file for each filter (and mask) in a chain. Therefore, it is good to always start by applying the `variable` filter where possible, and it may be necessary to clean unneeded files from `results/animal/vcf` occasionally.

###  Masks

The pipeline also has the concept of a mask. Whereas filters remove sites, masks set genotypes to missing. These include (curly braces here are [wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards)):

- `{dp}dp`: Mask genotypes with less than `{dp}` depth.
- `{het}het`: Mask genotypes with less than `{het}` heterozygous support.

These compose with filters in any order you like, e.g. `results/animal/vcf/Animal_A_variable_sites_10dp_3het_missing.bcf.gz`. Note that the order of the suffixes correspond to application order, from left to right. Hence, this would remove missing sites *after* masking depth and heterozygous support.

*Important*: Please note that masking may make various information in the masked variant record outdated. All masks update the `FORMAT/AC`, `FORMAT/AN`, and `FORMAT/MAF` tags after masking, but any other tags that are computed from genotype information will not be up-to-date. The [`fill-tags` plugin](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) for `bcftools` may be of use.

### Merging

When calling multiple groups separately, it is often required to merge the outgroups back into the ingroup. The pipeline treats this as somewhat akin to a left join on variants. In overview, all alleles observed in the ingroup will be kept, and any novel alleles introduced by the outgroups will be set to missing.

Note that unlike masks and filters, merging must take place after filtering variable sites and applying the sites filter. The target is then `results/animal/vcf/Animal_A_variable_sites_mergeB-C.bcf.gz` to merge outgroups `B` and `C` into the ingroup. See also the `merge_all_outgroups` pseudo-rule.

### Plink

The pipeline has limited support for conversion to [plink](https://www.cog-genomics.org/plink/1.9/formats) files. The naming follows the same rules as outlined above, but go in the `results/animal/plink` subdirectory and carry the usual plink extensions. To give an example, a valid target might be `results/animal/plink/Animal_A_variable_sites_mergeA-B.bed`.
