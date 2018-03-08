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
* [Command: overlap](#command-overlap)
  * [Overlap overview](#overlap-overview)
    * [Seed selection](#seed-selection)
    * [Indexing and querying](#indexing-and-querying)
    * [Chaining](#chaining)
  * [Overlap arguments](#overlap-arguments)
  * [Minimap2 performance comparison](#minimap2-performance-comparison)

# Installation
Downpore is written in Go and requires no external libraries. It has been tested on Linux but should compile and run on Mac OS and Windows too.

### Build from source
To build downpore from source:

* install [Go 1.8 or higher](https://golang.org/doc/install) and git
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
* `input` the input reads file
* `front_adapters` fasta containing adapters found at the beginning of reads
* `back_adapters` fasta containing adapters found at the end of reads
* `k` the size of k-mers to use in matching
* `chunk_size` the size to split long input reads into when indexing
* `middle_threshold` percentage identity to match an internal adapter
* `discard_middle` whether or not to keep the splits resulting from an internal match
* `check_reads` number of reads to consider when determining adapters present in the data
* `adapter_threshold` the identity to match an adapter when performing check_reads
* `extra_end_trim` number of bases to trim around adapters at the front and back
* `extra_middle_trim` number of bases to trim around internal adapters
* `tag_adapters` whether or not to modify output labels
* `verbosity` a level of logging output (written to stderr)
* `num_workers` number of threads to use to build the search index for internal adapters
* `himem` whether or not to cache input reads in memory

## De-multiplexing
The trim command does not explicitly de-multiplex data. If this is required then the sequence names of the output (which include barcode labels added by downpore) can be used by an external script to partition them.

## Porechop performance comparison
The main use case for the downpore trim command is for those situations in which Porechop is the bottleneck in your pipeline, or possibly when there are memory constraints. In terms of performance, downpore is I/O bound as it makes two passes through the input file, whereas Porechop is CPU bound. As such, all comparisons below are based on wall-clock time.

When the `--himem` argument is set the trim uses only a single pass through the input file and this has been tested separately.

### Performance test data
Three data sets: 
* E.coli small: 500MB of E.coli, R9.4 from Birmingham Uni (Bham_20171116_1xRAD004_6000ng.pass.fastq)
* E.coli: 1.5GB E.coli, R9.2 from Loman lab (E_coli_K12_1D_R9.2_SpotON_2.pass.fasta)
* Human: 2GB of chromosome 20 (na12878.chr20ScrappieFiltered.fasta) 

### Results

Trimming was performed on a low spec 4-core machine with an SSD.

Data | Speedup | Speedup (HiMem) | downpore memory | downpore (HiMem) memory | Porechop memory 
---| ---:| ---:| ---:| ---:| ---:
E.coli small | 31x | 34x | 0.3GB | 1.1GB | 1.1GB 
E.coli small .gz | 6x | 11x | 0.3GB | 1.1GB | 1.1GB 
E.coli | 25x | | 1.7GB | | 3.6GB 
E.coli .gz | 11x | | 1.7GB | | 3.6GB 
Human | 31x | | 1.3GB | | 4.6GB 
Human .gz | 12x | | 1.3GB | | 4.6GB 

In terms of adapters found, downpore typically finds a few percent more at the edges of reads. In the test examples, Porechop also applies a back adapter (which trimmed ~3-5% of reads) that is not clearly present in the first 10k reads but was included based on its association with a front adapter.

With downpore's `middle_threshold` set to 80 (rather than default 85) the two sets of splits are typically the same. Porechop reports a higher number of middle adapters as it includes some near the edges that downpore treats as a front/back trim instead.

For those repeating the above tests, also note that Porechop finds a number of false positive splits in E.coli. These come from an 86% identity match to the short back adapter that is actually present in the E.coli genome.

# Command: overlap
The overlap command is used to find overlaps amongst a set of long reads. This is similar to the function of [minimap2](https://github.com/lh3/minimap2) in `ava` mode, though the overlap command takes some shortcuts to improve performance.

Input reads are specified using the `-i` argument. Output overlaps are written to stdout in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).

Usage example:

```downpore overlap -i reads.fastq > overlaps.paf```

## Overlap overview
The overlap command matches the beginning and end of each input read (default length of 1000 bases) to the full set of input reads. Queries are performed in batches with each batch regenerating an index of all reads using a new set of seed k-mers. Sub-sequences with matching sets of k-mers then have these chained into gapped-seed sequences which are treated as a single ~1000 base overlap.

It is left to downstream analysis (assembly) to manage any ambiguous overlaps due to repeats in the genome. This includes the use of information pairing the left and right overlap results for each read. 

### Seed selection
By default each kilobase query has a minimum of 15 10-mer seeds. A batch of queries is complete once the limit of unique seeds is reached (default of 10000). Each query sequence generates new seeds until it reaches the required amount, so in practice early queries will have more matching seeds and the final query in the batch exactly 15 (or whatever the minimum is set to).

Seeds are chosen greedily so that they do not overlap in the query and minimise seed cost. The cost function is under development but currently:
* the lowest frequency (1%) of seeds are ignored
* the highest frequency (2%) are ignored
* seeds that appear in equal measure to their reverse complement are preferred

### Indexing and querying
The index contains a set of reads for each k-mer representing the "contains" relation. Long reads are split into overlapping chunks (default: 10000 base chunks) before being indexed. The index also contains the gapped-seed representation of each read.

When a query is made, the "soft-union" of all sets corresponding to k-mers present in the query is made. A soft-union with minimum `n` is the set of all elements that appear in at least `n` sets, e.g. the standard union operator is a soft-union with minimum 1. By default only reads that appear in at least 1/5th of all seeds' sets are sent on for chaining.

### Chaining
Every sequence returned from a query has its set of k-mers chained using dynamic programming. To be chained, k-mers must exactly match those of the query and the distance (in bases) between any adjacent k-mers must be within 66%-150% of their distance in the query.

## Overlap arguments
* `overlap_size` number of bases to use from the edge of each read
* `k` seed size in bases. Note: don't go above 13
* `num_seeds` maximum number of unique seeds in a query batch
* `seed_batch_size` maximum number of queries in a batch (if num_seeds not reached)
* `chunk_size` number of bases to split long reads into for indexing
* `num_workers` number of threads to use
* `input` the input reads file
* `himem` whether to cache input reads in memory (true by default)

## Minimap2 performance comparison
By default, `minimap2 -x ava-ont` uses the full input reads as queries. To make a fair comparison we extracted an query set that is equivalent to that used by downpore and ran `minimap2 -x ava-ont reads.fastq queries.fastq`. For small data with "normal" long reads this gives around a 3-4x speedup.

The minimap2 version used here was 2.7.

As ground truth we use minimap2 mapping to reference of both queries and all input reads. Generously, an overlap between a query read `q` and a target read `r` is considered correct if at least one (query) mapped location for `q` lies within at least one mapped region of `r`, and their complementarity matches.

At present I do not have a server for testing, so only the small E.coli data (as described in Command:trim) could be run under minimap2. The time and memory used for the standard minimap2 is listed as "ava" below in case it is of interest.

Algorithm | time | memory | recall | precision 
---| ---:| ---:| ---:| ---:| ---:
downpore | 35s | 0.7GB | 100%? | 95%?
minimap2 | 21s | 3.0GB | 83%?| 99.5%?
minimap2 (ava) | 1m31s | 3.1GB | | 

Downpore found strictly more overlaps than minimap2, however some of these were not identified by minimap2's map to reference. A manual inspection of a handful of these found some that appear to be correct overlaps -- so the 95% value listed is a lower bound on the true precision.

A true maximum for the recall value is not listed as the generous criteria used for counting correct overlaps would make this artificially high. The values shown are only useful as relative measures.
