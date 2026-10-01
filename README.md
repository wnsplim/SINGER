![Logo](SINGER.png)
# SINGER
SINGER stands for **S**ampling and **IN**ference of **GE**nealogies with **R**ecombination, and it is a Bayesian method to do posterior sampling of Ancestral Recombination Graph under Sequentially Markovian Coalescent. SINGER works by iterative threading one haplotype to the partially-built ARG, until the ARG for all haplotypes have been built. After initialization, MCMC will be performed to update the ARG to explore the posterior distribution. For a full description and cite our method, you can check: [Deng, Y., Nielsen, R., and Song, Y. S. "Robust and accurate Bayesian inference of genome-wide genealogies for hundreds of genomes." Nature Genetics 57, 2124–2135 (2025).](https://doi.org/10.1038/s41588-025-02317-9)

## Requirements

If you want to compile the source files, then C++17 and cmake are required. Otherwise you can also used the pre-compiled binary files on various platforms. 

The analysis on the inferred ARGs will be performed using [tskit](https://tskit.dev/tskit/docs/stable/introduction.html), and you can find the installation instructions using the link. **Note that now we support the compability with tskit version 1.0 onward, which requires python 3.10 and above (from version 0.1.9).** 

## Installations

The easiser way is to directory go to the folder `releases/` and download one of the versions which work for your working platform (Linux/MacOS_Intel/MacOS_M1). After downloading, you can decompress it using:

```
tar -xvzf file_name
```

To build from source, install zlib and htslib so that `pkg-config` finds htslib (add the directory that holds `htslib.pc` to `PKG_CONFIG_PATH` if needed), then run `bash beta_linux_compile.sh <version>` or `bash beta_mac_M1_compile.sh <version>` in `SINGER/SINGER`. The Python scripts read a `.bcf` input through `bcftools`, which must then be on `PATH`.

## Input and output

SINGER takes a **.vcf**, **.vcf.gz** or **.bcf** file and outputs a **.trees** file in tskit format. The mutations are already mapped to the branches, but non-polymorphic, multi-allelic sites and structral variants are excluded from inference. The branch length should be interpreted with units of generations, for example, for homo sapiens, you would need multiply that by 28 to convert to units of years. There will also be a **.log** file for you to check the argument you ran, and the summary statistic in MCMC iterations. 

## Basic usage

To sample ARGs with SINGER, you can run command line like shown below. 

**IMPORTANT!**:if you wish to get ARG for: **(1) a long chromosome or (2) a series of regions**, we have provided more support to help you (see the [next section: Tools](#Tools)). If you think there are other specific job pipeline which many people might want to use, please contact us and we might add it! 

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
|**-input**|required|input file: .vcf, .vcf.gz or .bcf|
|**-output**|required|prefix of the output .trees file name| 
|**-start**|required|start position of the region| 
|**-end**|required|end position of the region. SINGER reads the records with -start ≤ POS < -end.| 

The optional flags include:

|flag|required?|details|  
|-------------------|-----|---|  
|**-Ne**|optional|effective population size in diploid individuals: the number of haploid genomes is **2*Ne** at any ploidy, so a haploid population of N genomes takes -Ne N/2|
|**-ploidy**|optional|1 or 2, the ploidy of every sample in the VCF (default 2).|
|**-r**|optional|per base pair per generation recombination rate: a number, or the name of a rate map file. If `-m` is a number, the default is the value of `-m`. If `-m` is a rate map file, `-r` is required|
|**-n**|optional|the number of posterior samples, default at 100|
|**-thin**|optional|the number of MCMC iterations between adjacent samples, default at 200|
|**-polar**|optional|the probability of correct polarization, default at 0.5 for unpolarized data, please use 0.99 for polarized data|
|**-scaling_rep**|optional|the number of ARG rescaling rounds applied after the initial build and after each posterior sample, default at 1. Set to 0 to disable rescaling|
|**-scaling_bin**|optional|the number of time bins used for ARG rescaling, default at 100|
|**-tip_ages**|optional|file with the sampling ages in calendar years before present, see [Heterochronous samples](#heterochronous-samples-ancient-dna)|
|**-g**|optional|generation time in years, required with `-tip_ages`|
|**-mask**|optional|BED file of the regions that have no data, for example masked or unassayed regions. SINGER treats these regions as missing data and does not use the VCF records in them.|
|**-chrom**|optional|the chromosome to read from the input file (default: the first chromosome in the file)|

The output files will be:

```
prefix_of_output_files_nodes_{i}.txt, prefix_of_output_files_branches_{i}.txt, prefix_of_output_files_muts_{i}.txt, prefix_of_output_files_recombs_{i}.txt
```

with `i` from `0` to `num_samples - 1`. We recommend converting these files to tree sequence format in tskit, with this function:

```
path_to_singer/convert_to_tskit.py -input prefix_of_arg_files -output prefix_of_tskit_files
-start start_index -end end_index -step step_size
```

This tool will convert ARG sample with index from `start_index` to `end_index`, with interval size `step_size`. 

### Rate maps

`-m` and `-r` take a number or the name of a rate map file. A value that is a number is a constant rate. Any other value is the name of a rate map file. The two flags are independent, so you can give a recombination map with a constant mutation rate, or the opposite.

A rate map file has one segment on each line, with three numbers: the left position, the right position, and the rate per base pair per generation in that segment. The positions are chromosome positions, the same as in the VCF. Each segment must start where the previous segment ends. A rate of 0 is allowed. The map must cover the region from `-start` to `-end`. Example:

```
0 1000000 1.2e-8
1000000 2500000 0.8e-8
2500000 5000000 1.5e-8
```

With a recombination map, the bin width follows the local recombination rate, so recombination hotspots get narrower bins. Bins end at the segment boundaries of both maps.

SINGER stops with an error if the mutation rate is 0 everywhere between `-start` and `-end`, or if the VCF has a variant where the mutation rate is 0. If the recombination rate is 0 everywhere between `-start` and `-end`, SINGER prints a warning and infers an ARG without recombination.

With a mutation rate map, the ARG rescaling (`-scaling_rep`) weights the length of each branch by the mutation rate along its span.

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


## Tools

### Examining the convergence of the MCMC in SINGER

SINGER is an MCMC-based sampling algorithm. To examine the convergence of it we normally examine the traces of summary statistics, and we have found that 2 summary statistics are quite good at indicating the convergence of the SINGER MCMC: fit to the diversity landscape and number of non-uniquely-mapped sites, as used in the manuscript. We have provided a python script to calculate the traces for these 2 quantities:

```
python compute_traces.py
-prefix prefix_of_tree_sequence_file -m mutation_rate
-start_index index_start_sample -end_index index_terminal_sample
-output_filename output_trace_filename
```
The script will calculate the aforementioned two statistics for all samples with index between the `start_index` and the `end_index`, and they will be output as the the 2 columns in the output file. The fit to the diversity landscape and the number of non-uniquely-mapped sites typically drops with more iterations of the MCMC, and when they reach to a stable stage that is usually a good sign for convergence. 

### Computing the pairwise coalescence times with respect to a particular haplotype

In the manuscript we used the coalescence ratio to find introgression signals, which is based on the distribution of pairwise coalescence times between one haplytype and others, in every 10kb genome windows. Here we provide a python script tools, to compute the pairwise coalescence times between one given haplotype (indiciated by leaf node index) and others, in a windowed fashion. 

```
python compute_pairwise_coalescence_times.py
--trees_file tree_sequence_filename --leaf_index index_of_leaf_node
--interval_size size_of_genome_window --output_filename output_file_name
```
Each row in the output file stands for all the pairwise coalescence times between the input leaf node index and all others. The rows are in the order of the genome windows. 

### Running SINGER for a long chromosome

Often people would like to run the ARG inference method for the entire chromosome (or even the entire genome). We recommend running ```singer_master``` for continous segments (such as 5Mb) and then use the following tool to merge them together. 

```
path_to_singer/singer_master -m 1.25e-8 -input input.vcf.gz -output prefix_of_output_file_0 -start 0 -end 5e6
path_to_singer/singer_master -m 1.25e-8 -input input.vcf.gz -output prefix_of_output_file_1 -start 5e6 -end 10e6
path_to_singer/singer_master -m 1.25e-8 -input input.vcf.gz -output prefix_of_output_file_2 -start 10e6 -end 15e6
......
```

```
python merge_ARG.py --file_table sub_file_table_file --output merged_ARG_filename
```

The ```sub_file_table_file``` specifies how the inferred ARG should be pieced together with an example like this:

```
prefix_of_output_file_i_nodes_0.txt prefix_of_output_file_i_branches_0.txt prefix_of_output_file_i_muts_0.txt 0 
prefix_of_output_file_i_nodes_1.txt prefix_of_output_file_i_branches_1.txt prefix_of_output_file_i_muts_1.txt 5000000
prefix_of_output_file_i_nodes_2.txt prefix_of_output_file_i_branches_2.txt prefix_of_output_file_i_muts_2.txt 10000000
......
```

where ```i``` is the index of the ARG sample, and the index after that is the index of the genomic windows on which parallization is performed. Note that the genomic windows do not have to be all adjacent. For example, when there is a centromeric region in the middle of the chromosome, it can be skipped by only selecting windows in the left and right arms. This will be reflected in the last column of the ```sub_file_table_file```, which indicates the starting position of each genomic window. 

By doing the merging operation, you will get the ARG samples with the length of the entire chromosome rather than for each individual genomic window, which might be convenient in certain scenarios. 



## Suggestions from developer

1. As a Bayesian sampling method, SINGER works best when you sample some ARGs from posterior, **only using one single sample is NOT ideal**. To this point, we highly encourage specifying **-n, -thin** flags. You can find how we run SINGER on real datasets on:
2. It is of importance to carefully choose the parameters, such as -Ne, -m, and -r. We recommend first choosing the mutation rate m, and then based on average pairwise diversity \($\pi=4\cdot N_e \cdot m\$), you can decide the Ne parameter. If you are not super sure about the recombination rate, you can leave out -r, and SINGER sets it equal to -m. 
3. Unfortunately for now we only support phased, high-quality genomes, and polymorphic sites with missingness will be excluded. We are working on incorporating missingness and unphased data in the near future. ARGweaver has better support in these regards.
4. **A frequent bug reported** is caused by running SINGER on essentially a region with no or very low sequencing data, such as centromeric regions. SINGER cannot infer the ARG when there is no data present, and will likely bug out due to underflow issues. 
