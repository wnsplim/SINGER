![Logo](SINGER.png)
# SINGER
SINGER stands for **S**ampling and **IN**ference of **GE**nealogies with **R**ecombination, and it is a Bayesian method to do posterior sampling of Ancestral Recombination Graph under Sequentially Markovian Coalescent. SINGER works by iterative threading one haplotype to the partially-built ARG, until the ARG for all haplotypes have been built. After initialization, MCMC will be performed to update the ARG to explore the posterior distribution. For a full description and cite our method, you can check: [Deng, Y., Nielsen, R., and Song, Y. S. "Robust and accurate Bayesian inference of genome-wide genealogies for hundreds of genomes." Nature Genetics 57, 2124–2135 (2025).](https://doi.org/10.1038/s41588-025-02317-9)

## Requirements

If you want to compile the source files, then C++17 and htslib are required. Otherwise you can also use the pre-compiled binary files on various platforms. 

The analysis on the inferred ARGs can be performed using [tskit](https://tskit.dev/tskit/docs/stable/introduction.html), and you can find the installation instructions using the link. **Note that now we support the compatibility with tskit version 1.0 onward, which requires python 3.10 and above (from version 0.1.9).** 

The Python scripts need numpy, tskit, Cython and setuptools. `pyproject.toml` lists them and `uv.lock` pins their versions. To install them with [uv](https://docs.astral.sh/uv/), run `uv sync` in this folder, then run the scripts with `uv run python`.

## Installations

The easiest way is to directly go to the folder `releases/` and download one of the versions which work for your working platform (Linux/MacOS_Intel/MacOS_M1). After downloading, you can decompress it using:

```
tar -xvzf file_name
```

The Python scripts read a `.bcf` input through `bcftools`, which must then be on `PATH`.

## Input and output

SINGER takes a **.vcf**, **.vcf.gz** or **.bcf** file and outputs a **.trees** file in tskit format. Only single-nucleotide variants are used; other records (indels, structural variants, multi-allelic sites) are skipped, and regions covered by deletions are masked. The branch length should be interpreted with units of generations (unless explicit generation time is provided via the -g flag), for example, for *Homo sapiens*, you would need to multiply that by 28 to convert to units of years. There will also be a **.log** file with the summary statistics of the MCMC iterations.

## Basic usage

To sample ARGs with SINGER, you can run command line like shown below. 

**IMPORTANT!**:if you wish to get ARG for: **(1) a long chromosome or (2) a series of regions**, we have provided more support to help you (see [Running SINGER for a long chromosome](#running-singer-for-a-long-chromosome)). If you think there are other specific job pipeline which many people might want to use, please contact us and we might add it! 

```
path_to_singer/singer_master -m 1.25e-8
-input input.vcf.gz -output prefix_of_output_file
-start 0 -end 1e6
```

This command is to get the ARG samples for a specific region in the vcf file. We specify the details of the arguments here (or you can simply type ```path_to_singer/singer_master``` to display similar information):

The required flags include:

|flag|required?|details|  
|-------------------|-----|---|  
|**-m**|required|per base pair per generation mutation rate: a number, or the name of a rate map file (see [Rate maps](#rate-maps))|
|**-input**|required|input file: .vcf, .vcf.gz or .bcf. Without `-genotype_lik` or `-phase`, an unphased heterozygous genotype counts as missing for its sample.|
|**-output**|required|prefix of the output .trees file name| 
|**-start**|required|start position of the region| 
|**-end**|required|end position of the region. SINGER reads the records with -start ≤ POS < -end.| 

The optional flags include:

|flag|required?|details|  
|-------------------|-----|---|  
|**-Ne**|optional|effective population size: the diploid Ne under -ploidy 2, the haploid population size N under -ploidy 1. Without it `singer_master` estimates Ne from the diversity|
|**-ploidy**|optional|1 or 2, the ploidy of every sample in the input file (default 2).|
|**-r**|optional|per base pair per generation recombination rate: a number, or the name of a rate map file. If `-m` is a number, the default is the value of `-m`. If `-m` is a rate map file, `-r` is required|
|**-n**|optional|the number of posterior samples, default at 100|
|**-thin**|optional|the number of MCMC iterations between adjacent samples, default at 200|
|**-burnin**|optional|the number of MCMC iterations before the first sample, default at 2000. No sample written during the burn-in.|
|**-thin_log**|optional|the number of MCMC iterations between rows of the summary log, default at `-thin`/100. Set to 0 to write no summary log|
|**-thin_log_pos**|optional|the number of equidistant positions in the summary log, default at 100|
|**-polar**|optional|the probability of correct polarization, default at 0.5 for unpolarized data, please use 0.99 for polarized data|
|**-scaling_rep**|optional|the number of ARG rescaling rounds applied to a copy of each written posterior sample, default at 1; the rescaled node times are written to the sample file. Set to 0 to disable rescaling|
|**-scaling_bin**|optional|the number of time bins used for ARG rescaling, default at 100|
|**-mapping_rep**|optional|the number of stochastic mapping draws written for each posterior sample, default at 1. Each draw is written to the sample file. With 0 no draw is written.|
|**-tip_ages**|optional|file with the sampling ages in calendar years before present, see [Heterochronous samples](#heterochronous-samples-ancient-dna)|
|**-g**|optional|generation time in years, required with `-tip_ages`|
|**-mask_n**|optional|Negative mask. A BED file of the regions to mask; a row with a VCF sample name in a fourth column masks the region for that sample only.|
|**-mask_p**|optional|Positive mask. A BED file of the accessible regions, for example the 1000 Genomes strict accessibility mask for GRCh37 ([20141020.strict_mask.whole_genome.bed](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed)). You can use `-mask_p` together with `-mask_n`.|
|**-genotype_lik**|optional|use the GL or PL of each genotype, from the input record or else from `-genotype_calls`.|
|**-phase**|optional|sample the phase of each heterozygous genotype. Without it, a genotype written phased (`0\|1`) keeps its written phase; under `-genotype_lik`, an unphased heterozygous genotype has its phase sampled regardless of `-phase`.|
|**-genotype_calls**|optional|the caller's VCF whose GL or PL are used when the input has none, matched by position, REF, ALT and sample name. Read only with `-genotype_lik`.|
|**-chrom**|optional|the chromosome name to read from the input file (default: the first chromosome in the file)|
|**-chrom_n**|optional|the chromosome name to read from the `-mask_n` file (default: the `-chrom` chromosome)|
|**-chrom_p**|optional|the chromosome name to read from the `-mask_p` file (default: the `-chrom` chromosome)|
|**-resume**|optional|resume the interrupted run with the same `-output`. `-n` is the total number of posterior samples: the run continues from its last sample until `-n` samples exist. All other settings, including `-thin`, must match the original run.|

The output files will be:

```
prefix_of_output_files_{i}.bin
```

with `i` from `0` to `num_samples - 1`. The sample file holds the ARG of the sample, its rescaled node times when `-scaling_rep` is above 0, the counts, the mutations of each mapping draw, and the sampled genotypes under `-genotype_lik`. The counts are, for each edge: the expected number of mutations on the branch, and the number of its bases at which the data can see a mutation on it (a base where at least one sample below the branch and one sample outside it have data). Each mutation has: the position, the lower node of its branch, the state after it, and the informed flag (1 when the data can see it at that site by the same rule and 0 otherwise). When `-scaling_rep` is above 0, the mutations and counts are drawn at the rescaled node times. In Python, `read_sample` in `python_scripts/arg_io.py` reads the file, and `trees` yields the trees of the sample along the genome.

### Rate maps

`-m` and `-r` take a number or the name of a rate map file. The two flags are independent, so you can give a recombination map with a constant mutation rate, or the opposite.

A rate map file has one segment on each line, with three numbers: the left position, the right position, and the rate per base pair per generation in that segment. The positions are chromosome positions, the same as in the VCF. Each segment must start where the previous segment ends. A rate of 0 is allowed. The map must cover the region from `-start` to `-end`. Example:

```
0 1000000 1.2e-8
1000000 2500000 0.8e-8
2500000 5000000 1.5e-8
```


### Heterochronous samples (ancient DNA)

For data sets containing samples from different time points (e.g., ancient DNA), provide the sampling ages and generation time:

```
path_to_singer/singer_master -m mutation_rate -input input.vcf.gz -output prefix_of_output_file -start 0 -end 1e6 -tip_ages ages.txt -g 29
```

The tip ages file can be given in two formats, with the column count detected automatically. Ages are in calendar years before present, one per sample in the VCF (applied to both haplotypes of a diploid sample).

With **two columns** (`name  calendar_age_BP`), each age is matched by `name` to a sample name in the VCF header. Example:

```
tsk_0    0
tsk_1    3500
tsk_2    8000
```

With **one column**, give ages positionally in VCF sample order. The same three samples:

```
0
3500
8000
```

### Running SINGER for a long chromosome

Often people would like to run the ARG inference method for the entire chromosome (or even the entire genome). For this, `singer_master` cuts the region into blocks, runs SINGER on each block, and runs several blocks at the same time in parallel. Give the block length with `-L`:

```
path_to_singer/singer_master -m 1.25e-8 -input input.vcf.gz -output prefix_of_output_file
-start 0 -end 50e6 -L 5e6 -num_cores 10
```

Or give your own blocks in a file with `-blocks`:

```
path_to_singer/singer_master -m 1.25e-8 -input input.vcf.gz -output prefix_of_output_file
-start 0 -end 50e6 -blocks blocks.txt -num_cores 10
```

Each line of the block file is one block: its left and its right position, in the chromosome positions of the VCF, as in a rate map. The lines are sorted, the blocks do not overlap, and all blocks are between `-start` and `-end`. Gaps between the blocks are allowed. For example, these blocks leave out a region from 18 Mb to 22 Mb:

```
0 6000000
6000000 12000000
12000000 18000000
22000000 28000000
```

|flag|details|
|-------------------|---|
|**-L**|block length: `singer_master` cuts the region from `-start` to `-end` into blocks of this many bases|
|**-blocks**|file of blocks, one `left right` line for each block|
|**-num_cores**|the number of blocks that run at the same time, default at 1|

Block `b` (counted from 0, in the order of the blocks) writes the output files of a SINGER run with the prefix `prefix_of_output_file_block{b}`, and its console output to `prefix_of_output_file_block{b}.out`. It runs with the seed `-seed` + `b`. In a file name such as `prefix_of_output_file_block{b}_{i}.bin`, `b` is the block index and `i` is the sample index; the converted `prefix_of_output_file_{i}.trees` is sample `i` of the whole region. A block that holds no VCF record is skipped. Without `-Ne`, `singer_master` estimates Ne one time from the diversity of all the blocks that run, with the masks applied, and gives this value to every block. All the other flags apply to every block.

The file `prefix_of_output_file.blocks` lists the region and the blocks that ran. The first line has the start and the end of the region. Each other line is one block: its index, its left position and its right position. For the first command above, when the block from 15 Mb to 20 Mb holds no VCF record:

```
0.0	50000000.0
0	0.0	5000000.0
1	5000000.0	10000000.0
2	10000000.0	15000000.0
4	20000000.0	25000000.0
......
```

### Converting to tskit

We recommend converting these files to tree sequence format in tskit, with this function:

```
path_to_singer/convert_to_tskit.py -input prefix_of_output_file.ckp -output prefix_of_tskit_files
-start start_index -end end_index -step step_size
```

|flag|details|
|-------------------|---|
|**-input**|the `.ckp` file of the SINGER run, or its `.blocks` file for a run in blocks|
|**-output**|prefix of the .trees files to write|
|**-start**|index of the first posterior sample to convert|
|**-end**|index "after" the last posterior sample to convert|
|**-step**|spacing between the converted sample indices, default at 1; for example `-start 0 -end 100 -step 10` converts samples 0, 10, ..., 90|
|**-vcf**|optional: the input VCF of the run, whose sample names are attached to the tree sequence as individuals|

Sample `i` is written to `prefix_of_tskit_files_{i}.trees`. With `-mapping_rep` above 1 the tool writes one tree sequence for each draw, `prefix_of_tskit_files_{r}_{i}.trees`. The tool prints the name of each file it writes.

Each tree sequence holds:

- the ARG of the sample, with the rescaled node times when SINGER wrote them, else the sampled ones;
- the mutations of the mapping draw, with the `informed` flag in the mutation metadata; their times are computed by tskit (`compute_mutation_times`). Without a mapping draw (`-mapping_rep 0`) the tree sequence will not carry any mutation table;
- the counts in its edge metadata (`expected_mutations`, `informed_bases`). Call `tables.edges.drop_metadata()` before `simplify`, which tskit refuses on edges with metadata;
- the flags of the run in its provenance record.

For a run in blocks (`-L` or `-blocks`), the same command with `-input prefix_of_output_file.blocks` stitches the blocks into one tree sequence. `convert_to_tskit.py` places each block at its chromosome position in one tree sequence. The sequence length is the end of the region. A gap between the blocks, and a block that was skipped, has no edges and no sites. To leave a block out of the tree sequence, delete its line from `prefix_of_output_file.blocks` before you convert.

You can also give `-convert_step k` to `singer_master`. At the end of the run, it then converts every `k`-th sample to `prefix_of_output_file_{i}.trees` with the same content, and attaches the sample names of the input file.


## Tools

### Examining the convergence of the MCMC in SINGER

SINGER is an MCMC-based sampling algorithm. To examine the convergence of it, SINGER writes the summary log `prefix_of_output_files_summary.log`: one row every `-thin_log` MCMC iterations, from the start of the burn-in, with the number of recombinations, and at each of the `-thin_log_pos` equidistant positions the total branch length and the 25, 50 and 75% quantiles of the node ages of the local tree. The local trees at these positions are written to `prefix_of_output_files_summary_trees.bin`. We have provided a python script to compute the effective sample sizes (ESS) from these two files:

```
python compute_ess.py -input prefix_of_output_files_summary.log -output prefix_of_ess_file
[-burnin burn_in_iterations] [-tree_ess rf rnni triplet]
```

|flag|details|
|-------------------|---|
|**-input**|the summary log of the SINGER run; the script reads the trees file and the `.ckp` file of the run from the same folder|
|**-output**|prefix of the .tsv file to write|
|**-burnin**|the rows of the summary log up to this MCMC iteration are left out, default at the `-burnin` of the run|
|**-tree_ess**|the tree distances for the tree ESS, one or more of `rf` (Robinson-Foulds), `rnni` (ranked nearest neighbour interchange) and `triplet`, default at `rf`; giving `0` skips all tree ESS computation|

The script computes the ESS of every column of the summary log, and for the local tree at each position the minimum and median pseudo-ESS and the Fréchet correlation ESS of [Magee et al. (2024)](https://doi.org/10.1214/22-BA1339). It prints their quantiles over the positions and writes every value to `prefix_of_ess_file_ess.tsv`, one row per column of the summary log, then one row per tree ESS named `{measure}_{distance}_{position index}`. The tree ESS needs Cython and a C compiler, and the script compiles `tree_distances.pyx` on its first run.

### Computing the pairwise coalescence times with respect to a particular haplotype

In the manuscript we used the coalescence ratio to find introgression signals, which is based on the distribution of pairwise coalescence times between one haplotype and others, in every 10kb genome windows. Here we provide a python script tools, to compute the pairwise coalescence times between one given haplotype (indicated by leaf node index) and others, in a windowed fashion. 

```
python compute_pairwise_coalescence_times.py
-input tree_sequence_filename -leaf_index index_of_leaf_node
-interval_size size_of_genome_window -output prefix_of_output_file
```
Each row in the output file `prefix_of_output_file_pairwise_coalescence_times.csv` stands for all the pairwise coalescence times between the input leaf node index and all others. The rows are in the order of the genome windows. 

## Suggestions from developer

1. As a Bayesian sampling method, SINGER works best when you sample some ARGs from posterior, **only using one single sample is NOT ideal**. To this point, we highly encourage specifying **-n, -thin** flags. You can find how we run SINGER on real datasets on:
2. It is of importance to carefully choose the parameters, such as -Ne, -m, and -r. We recommend first choosing the mutation rate m, and then based on average pairwise diversity \($\pi=4\cdot N_e \cdot m\$ under -ploidy 2, $\pi=2\cdot N \cdot m\$ under -ploidy 1), you can decide the Ne parameter. If you are not super sure about the recombination rate, you can leave out -r, and SINGER sets it equal to -m. 
3. **A frequent bug reported** is caused by running SINGER on essentially a region with no or very low sequencing data, such as centromeric regions. SINGER cannot infer the ARG when there is no data present, and will likely bug out due to underflow issues. 
