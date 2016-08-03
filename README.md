# damapper_bwt
damapper plus BWT based seed searching

Source
------

The damapper_bwt source code is hosted on github:

        https://github.com/gt1/damapper_bwt

It can be obtained using

```
git clone --recursive https://github.com/gt1/damapper_bwt
```

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

The program expects the reference and reads in the FastA format [https://en.wikipedia.org/wiki/FASTA_format]
and outputs alignments as an (uncompressed) BAM file on the standard output channel.

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
