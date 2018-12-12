# downpore
A suite of tools for use in genome assembly and consensus. Work in progress.

# Table of contents
* [Installation](#installation)
  * [Build from source](#build-from-source)
  * [Run precompiled binary](#run-precompile-binary)
* [Usage](#usage)
* [Command: trim](#command-trim)
  * [Trim overview](#trim-overview)
    * [Matching criteria](#matching-criteria)
    * [Sequence labels](#sequence-labels)
  * [Trim arguments](#trim-arguments)
  * [De-multiplexing](#de-multiplexing)
  * [Porechop performance comparison](#porechop-performance-comparison)
* [Command: map](#command-map)
  * [Map overview](#map-overview)
    * [Seed selection for mapping](#seed-selection-for-mapping)
    * [Indexing and querying](#indexing-and-querying)
    * [Chaining](#chaining)
  * [Map arguments](#overlap-arguments)
  * [Minimap2 mapping comparison](#minimap2-mapping-comparison)
    * [Mapping precision and recall](#mapping-precision-and-recall)
    * [Mapping run time](#mapping-run-time)
* [Command: overlap](#command-overlap)
  * [Overlap overview](#overlap-overview)
    * [Seed selection for overlaps](#seed-selection-for-overlaps)
  * [Overlap arguments](#overlap-arguments)
  * [Minimap2 overlap comparison](#minimap2-overlap-comparison)
* [Command: subseq](#command-subseq)
* [Evaluation data](#evaluation-data)

# Installation
Downpore is written in Go and requires no external libraries. It has been tested on Linux but should compile and run on Mac OS and Windows too.

### Build from source
To build downpore from source:

* install [Go 1.9 or higher](https://golang.org/doc/install) and git
* create a directory `$GOPATH/src/github.com/jteutenberg`
* from that directory, clone the downpore repository `git clone https://github.com/jteutenberg/downpore.git`
* enter the new downpore directory and build it using `go build downpore.go`

### Run precompiled binary
A stand-alone, pre-compiled binary for Linux x86-64 is available at https://github.com/jteutenberg/downpore/releases/latest

# Usage
The general form of execution is
```downpore <command> [arguments]```
where arguments are provided with either a one- or two-hyphen switch followed by a space or equals sign and then its value.

Running `downpore` with no arguments gives a list of available commands.

To see the available arguments for a command, use `downpore help <command>`.

## Read input format
Input reads must be a fasta/fastq file or a gzip (with .gz suffix) containing a fasta/fastq. Fasta files with multiple lines per sequence are not handled and will need to have end-lines removed.

Only DNA sequences composed of A,C,G and T are accepted. Other characters will be assigned one of these values.

# Command: trim
The trim command is used to remove adapters or barcodes from the end of long reads. It broadly performs the same function as [Porechop](https://github.com/rrwick/Porechop), and while downpore is less polished it is substantially faster.

Input reads are specified using the `-i` argument. Output reads are written to stdout and will be in the same format as the input reads.

Usage examples:

```downpore trim -i reads.fastq -f ./data/adapters_front.fasta -b ./data/adapters_back.fasta > trimmed.fastq```

```downpore trim -i reads.fastq -f ./data/adapters_front.fasta -b ./data/adapters_back.fasta --himem true --num_workers 32 > trimmed.fastq```

```downpore trim -i reads.fastq -f ./data/adapters_front.fasta -b ./data/adapters_back.fasta --num_workers 32 --demultiplex ./barcode_output --require_pairs true```

## Trim overview
The trim command uses k-mer matching and chaining to find sub-sequences of reads that match any adapter/barcode from a list provided by the user.  Adapters found in the middle of long reads cause that read to be split.

Example adapter lists can be found in `data/adapters_front.fasta` and `data/adapters_back.fasta`.

Where possible, arguments and their default values mirror those of Porechop. Important functional differences are that downpore trim does not pair up front and back adapters (all are treated as being independent); and there is no notion of good and bad sides of an internal adapter: both are conservatively treated as bad.

### Matching criteria
Adapters are identified with high recall at the edges (first and last 150 bases) and with higher precision in the middle of a read.

By default, 6-mers are used in matching. At the edges, an adapter is identified as being present when at least 3 k-mers match in-order and at approximately the expected distance from one another. This means that at a minimum 8 contiguous matching bases are present.

Internal matches are made when the "identity" matching passes the threshold specified by `-middle_threshold` (default 85%). The identity value is taken as the percentage of bases in the adapter that are contained in at least one matching k-mer. 

### Sequence labels
By default, the adapter/barcode with the most bases present in a read has its name appended to the beginning of the read's label in the trimmed output. This can be turned off using `-tag_adapters`.

Adapters with names beginning "Barcode" are a special case. These take precedence and will always be used in the output label if found, though the trimming will still be based on all adapters present. If there exist two barcodes that are within 5% identity in a read (i.e. ambiguous barcodes) then no adapter labels will be written. Take care: this has not been thoroughly tested.

## Trim arguments
Input
* `input` the input reads file
* `front_adapters` fasta containing adapters found at the beginning of reads
* `back_adapters` fasta containing adapters found at the end of reads

Adapter sub-sets to use
* `determine_adapters` whether to use all adapters provided or to determine which are present
* `check_reads` number of reads to consider when determining adapters present in the data
* `adapter_threshold` the identity to match an adapter when performing check_reads

Matching parameters
* `k` the size of k-mers to use in matching
* `chunk_size` the size to split long input reads into when indexing
* `middle_threshold` percentage identity to match an internal adapter
* `discard_middle` whether or not to keep the splits resulting from an internal match
Trim output parameters
* `extra_end_trim` number of bases to trim around adapters at the front and back
* `extra_middle_trim` number of bases to trim around internal adapters
* `tag_adapters` whether or not to modify output labels

Additional barcode parameters
* `demultiplex` a path to demultiplex barcodes into, otherwise write everything to stdout
* `require_pairs` whether to ignore reads without a matching pair of barcodes at both ends

Other arguments
* `verbosity` a level of logging output (written to stderr)
* `num_workers` number of threads to use to build the search index for internal adapters
* `himem` whether or not to cache input reads in memory

## De-multiplexing
When the demultiplex path is provided sequences will be written to one fasta/fastq file per barcode. A barcode is any adapter that begins with "Barcode" and should not contain any underscores in its name. The `tag_adapters` switch should be left true, though in this case the demultiplexed sequences will not have a prefix added to their names.

When a large number of barcodes is in use the set of k-mers covered by the full set will become large. This can greatly increase memory usage. Should this become an issue the simplest solution is to increase k (say, to k=7 from the default of 6) so that their presence in the reads is less dense.

## Porechop performance comparison
The main use case for the downpore trim command is for those situations in which Porechop is the bottleneck in your pipeline, or possibly when there are memory constraints. In terms of performance, downpore is I/O bound as it makes two passes through the input file, whereas Porechop is CPU bound. As such, all comparisons below are based on wall-clock time.

When the `--himem` argument is set the trim uses only a single pass through the input file and this has been tested separately.

See the [evaluation data](#evaluation-data) section for details on the computer and data sets used below.

Data | Speedup | Speedup (HiMem) | downpore memory | downpore (HiMem) memory | Porechop memory 
---| ---:| ---:| ---:| ---:| ---:
E.coli small | 31x | 34x | 0.3GB | 1.1GB | 1.1GB 
E.coli small .gz | 6x | 11x | 0.3GB | 1.1GB | 1.1GB 
E.coli | 25x | | 1.7GB | | 3.6GB 
E.coli .gz | 11x | | 1.7GB | | 3.6GB 
Human Ch20 | 31x | | 1.3GB | | 4.6GB 
Human Ch20 .gz | 12x | | 1.3GB | | 4.6GB 

In terms of adapters found, downpore typically finds a few percent more at the edges of reads. In the test examples, Porechop also applies a back adapter (which trimmed ~3-5% of reads) that is not clearly present in the first 10k reads but was included based on its association with a front adapter.

With downpore's `middle_threshold` set to 80 (rather than default 85) the two sets of splits are typically the same. Porechop reports a higher number of middle adapters as it includes some near the edges that downpore treats as a front/back trim instead.

For those repeating the above tests, also note that Porechop finds a number of false positive splits in E.coli. These come from an 86% identity match to the short back adapter that is actually present in the E.coli genome.

# Command: Map
The map command is used to find the most likely approximate source location in a reference for each input read. This is analogous to approximate long read mappers like [minimap2](https://github.com/lh3/minimap2).

Input reads are specified using the `-i` argument and a single target reference sequence as a fasta using the `-r` argument. Output mappings are written to stdout in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).

Usage example:

```downpore map -i reads.fastq -r reference.fasta -circular true > mappings.paf```

## Map overview
The map command first creates an index of the reference then queries portions of each input read (default 1k bases) against it, working in from the edges. Both index and queries are represented by a sequence of seed k-mers (default 11-mers) and the size of the gaps between them. Initially, the set of seeds in a query are matched against the index, then these are chained to find the best in-order matches.

There is no restriction on the number of matches output, though clearly lower identity matches that are subsequences of other results are discarded.

During querying, when a pair of consistent matches is found for both ends of a query sequence a candidate output match is generated. If no such pairs are found the query is assumed to be chimeric and further mappings in the form of a binary search are perfomed to find the split point.

### Seed selection for mapping
Seeds are selected from the reference sequence, with a minimum of one seed per 75 bases by default. For each 75-base subsequence, if no prior selected seeds are present, a new seed is chosen to maximise the following criteria. This is a work in progress:
* the lowest frequency (1%) of seeds are ignored
* the highest frequency (2%) are ignored
* seeds with frequency approaching 1 per 500k are preferred

The preferred seed frequency is based on seeds of 10-12 bases. A balance needs to be struck between rare seeds (good for reducing false positives) and common seeds (improving efficiency by reducing the total seed count).

### Indexing and querying
The index has a set of reads for each k-mer representing the "contains" relation. Long reads are split into overlapping chunks (default: 10000 base chunks) before being indexed. The index also contains the gapped-seed representation of each read.

When a query is made, the "soft-union" of all sets corresponding to k-mers present in the query is made. A soft-union with minimum `n` is the set of all elements that appear in at least `n` sets, e.g. the standard union operator is a soft-union with minimum 1. By default only reads that appear in at least 1/4 of all seeds' sets are sent on for chaining.

### Chaining
Every sequence returned from a query has its set of k-mers chained using dynamic programming. To be chained, k-mers must exactly match those of the query and the distance (in bases) between any adjacent k-mers must be within 66%-150% of their distance in the query.

## Map arguments
* `k` seed size in bases. Note: don't go above 13
* `query_size` number of bases to use in each sub-query
* `chunk_size` number of bases to split the reference into for indexing
* `seed_rate` minimum number of bases per seed
* `num_workers` number of threads to use
* `input` the input reads file
* `reference` the reference sequence file
* `circular` whether this reference is for a circular genome

## Minimap2 mapping comparison
Here we compare downpore's map with minimap2 (version 2.7) run as `minimap2 -x map-ont -N 1000` which should have comparable behaviour.

The "ground truth" is taken as the union of mappings made by the two mappers, with the mapping identity calculated by Needleman-Wunsch using the `needle` command from EMBOSS. The large E.Coli and human Ch20 datasets described in the [evaluation data](#evaluation-data) section are used here.

### Mapping precision and recall

To begin we consider the overlap between the algorithms: the number of mappings shared by both, those found only by downpore, and found only by minimap2. Here only mappings to more than 500 bases of the reference are included.

Dataset | found by both | downpore only | minimap2 only
---| ---:| ---:| ---:
E.Coli | 154k | 5k | 3k 
Human Ch20 | 245k | 397k | 129k

Note that the human chromosome 20 contains repeats. For example, downpore found unique mapping locations for 219k sequences, and multiple mapping locations for only 21k sequences, but each of the multiple-location hits averaged around 20 mappings per sequence.

The mappings found above are the combination of both true and false positives. We take a "false positive" to be either a mapping with very low identity (say, below 50%) or with identity substantially lower (more than 15% lower) than the best identity for all mappings matched using that subsequence. The two plots below show the distribution of reads' identities for the three categories shown in the table above, followed by the distribution of differences from the best identities for the same set of reads.

<img src="https://github.com/jteutenberg/downpore/blob/master/img/docs_mapping_ecoli_identity.png" width="400"><img src="https://github.com/jteutenberg/downpore/blob/master/img/docs_mapping_human_identity.png" width="400">

On E.coli minimap2 finds a few additional mappings that are slightly higher quality than those found by downpore. Conversely, on the repetitive human genome downpore ignores the low identity population of mappings while finding more overall.

<img src="https://github.com/jteutenberg/downpore/blob/master/img/docs_mapping_ecoli_deltaidentity.png" width="400"><img src="https://github.com/jteutenberg/downpore/blob/master/img/docs_mapping_human_deltaidentity.png" width="400">

The relative identities shown in the distributions above can be used to determine the precision by assigning false/true positives. The *mappings* recalled by at least one mapper are used to determine the "total recall". Later we will also calculate a recall value over the *input sequences* rather than the mappings (i.e. whether at least one valid mapping was found) and we use the term total recall here to distinguish from this.

Dataset | downpore total recall | minimap2 total recall | downpore precision | minimap2 precision
---| ---:| ---:| ---:| ---:
E.Coli | 98.6% | 97.8% | 99.98% | 99.98% 
Human Ch20 | 79.9% | 40.4% | 99.2% | 96.9%

On E.coli the algorithms have similar results. On the human ch20 data downpore has greater recall and the additional results it finds tend to be slightly higher identity than the additional mappings found by minimap2. 

One fact not shown in the plots is that 19k of the mappings found by both algorithms have a longer version found by minimap 2 -- i.e. the reference sequence identified by downpore is a subsequence of that found by minimap2. Conversely, only 5k mappings from minimap2 were subsequences of those of downpore. So it can be said that minimap2 is better able to extend mappings out.

Finally, the low recall of minimap2 is due to the high number of repeat matches. Some use cases are only interested in the presence of at least one good match in which case the values above are not relevant. Below are the (non-total) recall values:

Dataset | downpore recall | minimap2 recall
---| ---:| ---:|
E.Coli | 99.9% | 99.7%
Human Ch20 | 98.4% | 97.0%

With both mappers finding the best mapping for 95.4% of the human sequences, and 99.5% of the E.coli. So for finding the top match these are more or less equivalent. 

### Mapping run time

Dataset | downpore time | minimap2 time | downpore memory | minimap2 memory |
---| ---:| ---:| ---:| ---:|
E.Coli | 6.7s | 10.7s | 0.4GB | 1.0GB
Human Ch20 | 48.7s | 33.5s | 1.0GB | 1.9GB

Downpore is more parsimonious with memory. In the presence of repeats it takes longer than minimap2, with mostly unique mappings it is faster on these fairly small sets of data.

# Command: overlap
The overlap command is used to find overlaps amongst a set of long reads. This is similar to the function of [minimap2](https://github.com/lh3/minimap2) in `ava` mode, though the overlap command takes some shortcuts to improve performance.

Input reads are specified using the `-i` argument. Output overlaps are written to stdout in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).

Usage example:

```downpore overlap -i reads.fastq > overlaps.paf```

## Overlap overview
The overlap command matches the beginning and end of each input read (default length of 1000 bases) to the full set of input reads. Queries are performed in batches with each batch regenerating an index of all reads using a new set of seed k-mers. Sub-sequences with matching sets of k-mers then have these chained into gapped-seed sequences which are treated as a single ~1000 base overlap.

It is left to downstream analysis (assembly) to manage any ambiguous overlaps due to repeats in the genome. This includes the use of information pairing the left and right overlap results for each read. 

### Seed selection for overlaps
By default each kilobase query has a minimum of 15 10-mer seeds. A batch of queries is complete once the limit of unique seeds is reached (default of 10000). Each query sequence generates new seeds until it reaches the required amount, so in practice early queries will have more matching seeds and the final query in the batch exactly 15 (or whatever the minimum is set to).

Seeds are chosen greedily so that they do not overlap in the query and minimise seed cost. The cost function is under development but currently:
* the lowest frequency (1%) of seeds are ignored
* the highest frequency (2%) are ignored
* seeds that appear in equal measure to their reverse complement are preferred

## Overlap arguments
* `overlap_size` number of bases to use from the edge of each read
* `k` seed size in bases. Note: don't go above 13
* `num_seeds` maximum number of unique seeds in a query batch
* `seed_batch_size` maximum number of queries in a batch (if num_seeds not reached)
* `chunk_size` number of bases to split long reads into for indexing
* `num_workers` number of threads to use
* `input` the input reads file
* `himem` whether to cache input reads in memory (true by default)

## Minimap2 overlap comparison
By default, `minimap2 -x ava-ont` uses the full input reads as queries. To make a fair comparison we extracted an query set that is equivalent to that used by downpore and ran `minimap2 -x ava-ont reads.fastq queries.fastq`. For small data with "normal" long reads this gives around a 3-4x speedup.

An overlap between a query read `q` and a target read `r` is considered correct if at least one (query) mapped location for `q` lies within at least one mapped region of `r` (according to the "ground truth" determined in the map command section above), and their complementarity matches.

The time and memory used for the standard minimap2 is also listed as "ava" below in case it is of interest.

### E.Coli

These tables will be regenerated soon using the ground truth used in the map command section.

Algorithm | time | memory | recall | precision 
---| ---:| ---:| ---:| ---:| ---:
downpore | | | | 
minimap2 | | | |
minimap2 (ava) | | | | 

### Human Ch20

Algorithm | time | memory | recall | precision 
---| ---:| ---:| ---:| ---:| ---:
downpore | | | | 
minimap2 | | | |
minimap2 (ava) | | | | 

# Command: subseq
The subseq commands starts a process that, given an input fasta, accepts requests from stdin for subsequences and outputs the result to stdout. In the first pass through the file the offsets of each sequence within the file are stored (or the entire sequence stored in memory) to enable fast retrieval.

Usage example:

```downpore subseq -i sequences.fasta```

Commands to stdin are of the form:

```<start base> <end base> <true/false> [sequence name]```

where the true/false specify whether or not a reverse-complement output is required, and the sequence name specifies the entry in the input fasta file. If no name is given then the subsequence is retrieved from the first sequence in the input file.

Output is either two lines in fasta format, or a single line error message.

# Evaluation data

Three data sets are used for evaluating the commands listed above: 
* E.coli small: 500MB of E.coli, R9.4 from Birmingham Uni (Bham_20171116_1xRAD004_6000ng.pass.fastq)
* E.coli: 1.5GB E.coli, R9.2 from Loman lab (E_coli_K12_1D_R9.2_SpotON_2.pass.fasta)
* Human: 2GB of chromosome 20 (na12878.chr20ScrappieFiltered.fasta) 

All evaluations are run on a single 8-core/16-thread machine, with programs set to use 16 threads when this option is available. Data is stored on an M.2 SSD and memory clocked at 2400 MHz.

All times are measured as wall-clock and memory measurements are the peak usage throughout execution.
