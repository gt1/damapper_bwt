# damapper_bwt
damapper plus BWT based seed searching

Source
------

The damapper_bwt source code is hosted on github:

        https://github.com/gt1/damapper_bwt

The current development branch can be obtained via

```
git clone --recursive https://github.com/gt1/damapper_bwt
```

Release tar balls can be found at

	https://github.com/gt1/damapper_bwt/releases

Please check for versions with -release- in their name.

Binary releases
---------------

Binary releases for Linux x86_64 can be found at https://github.com/gt1/damapper_bwt/releases/ .

Compilation of damapper_bwt
---------------------------

damapper_bwt needs libmaus2 [https://github.com/gt1/libmaus2].
When libmaus2 is installed in ${LIBMAUSPREFIX} then damapper_bwt can be compiled and
installed in ${HOME}/damapper_bwt using

        - autoreconf -i -f
        - ./configure --with-libmaus2=${LIBMAUSPREFIX} --prefix=${HOME}/damapper_bwt
        - make install

Running damapper_bwt
--------------------

damapper_bwt can be called using

```
damapper_bwt ref.fasta <reads.fasta >mapped.bam
```

The program expects the reference in FastA format [https://en.wikipedia.org/wiki/FASTA_format].
The reads can be provided in FastA or BAM format. If the input is in the BAM format then the program needs to be called using -Ibam.
damapper_bwt outputs alignments as a BAM file on the standard output channel.

During the first run with any reference the program will compute an index of the reference.

The following parameters can be used (no space allowed between option an argument):

* -k: seed length (default 20)
* -K: kmer cache word length (default 12)
* -p: number of threads (defaults to number of cores on machine)
* -v: verbosity level (default: 1)
* -i: maximum number of read bases per input block (default 64m)
* -M: damapper memory limit
* --sasamplingrate: SA sampling rate (default 32)
* --bwtconstrmem: memory used to construct BWT (default 3/4 of machine's memory)
* -T: prefix for temporary files during index construction
* -Q: index file name (defaults to .damapper_bwt appended to reference file name)
* -S: sequence storage strategy for non primary alignments (none (default), soft, hard)
* -z: output BAM compression level (zlib default)
* -I: input format (fasta (default) or bam)

BAM output files
----------------

damapper_bwt outputs alignments as a BAM file. As damapper deliberately (see https://dazzlerblog.wordpress.com/2016/07/31/damapper-mapping-your-reads/)
outputs more than one alignment per read for several reasons like

1. there is more than one possible locus to place a read on the reference
2. a read contains low quality regions which do not align well, so damapper chains up the good stretches while leaving the bad ones unaligned

this will be reflected in the BAM file as well. The alignments of the highest scoring chain output by damapper for a read do not have the secondary flag set in the BAM file.
All but the first alignment of a damapper chain have the supplementary flag set in the BAM file. The following auxiliary fields can be used to deduce the chain id of an alignment and
the position of an alignment inside its chain:

* cn (integer) stores the number of chains reported for the query sequence/read
* ci (integer) stores the chain id an alignment belongs to (this ranges from 0 to the value stored in cn minus one)
* cl (integer) stores the length of the chain this alignment belongs to (i.e. the chain is comprised of this many alignments)
* cj (integer) stores the chain link id of the alignment (this ranges from 0 to the value stored in cl minus one)

Note that only the first alignment line for a read (the one with ci=0 and cj=0 if the read is aligned) store the query sequence when the program is
run without using the -S switch or if the switch is set to -Snone. All other lines for a read have an undefined query sequence field.
Some downstream programs cannot handle alignments which do not store the query sequence. This should be taken into account when
calling damapper_bwt. The best compatibility (but also the biggest output files) is obtained by using -Ssoft. In this mode
each output alignment contains the complete query sequence.

If the input is given in BAM format, then damapper_bwt will copy read group (RG) identifiers and other auxiliary tags from
the input to the output file.
