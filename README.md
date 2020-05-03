## GffCompare      

* compare and evaluate the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). 
* collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples)
* classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a
annotation file (also in GTF/GFF3 format)

More details and usage examples can be found in the paper: [DOI: 10.12688/f1000research.23297.1](http://dx.doi.org/10.12688/f1000research.23297.1) which can be also used to cite this software.

The official documentation and download packages for this utility can be found online here:
http://ccb.jhu.edu/software/stringtie/gffcompare.shtml

For more information about the GFF3/GTF file formats expected by
GffCompare please refer to the available online documentation -- a quick
review can be found at: http://ccb.jhu.edu/software/stringtie/gff.shtml

The original version of this program was distributed as part of the
Cufflinks suite, under the name "CuffCompare".

The overall functionality and most of the options of CuffCompare are
still supported by GffCompare, while new functionality is being added to
GffCompare only, as it is the program which is actively maintained.

An example of a new feature of GffCompare (compared to its predecessor
CuffCompare) is this: when a single query GTF/GFF file is given as input
for analysis, along with a reference annotation (-r option), GffCompare
switches into *annotation mode* and it generates a *.annotated.gtf* file
instead of the *.combined.gtf* produced by CuffCompare with the same
parameters. This file has the same general format as CuffCompare's
*.combined.gtf* file (with "class codes" assigned to transcripts as per
their relationship with the matching/overlapping reference transcript),
but the original transcript IDs are preserved, so GffCompare can thus be
used as a simple way of annotating a set of transcripts.

Another important difference is that the input transcripts are by default no longer discarded when they are found to be "intron redundant", i.e. contained within other, longer isoforms. CuffCompare had the -G option to prevent collapsing of such intron redundant isoforms into their longer "containers", but GffCompare has made this the default mode of operation (hence the -G option is no longer needed and is simply ignored when given). However please note that "matching" transcripts with fully identical intron chains (i.e. with the same exact intron coordinates, hence the same intron-exon structure except the terminal exon ends) are *still* discarded when GTF/GFF files are loaded.

## trmap
Some pipelines can produce a very large number of potential or partial transcripts ("transfrags"), for example when merging the transcript assemblies from tens or hundreds of RNA-Seq experiments assemblies with `stringtie --merge`. Running GffCompare on such large GTF/GFF files could be slow and memory intensive (because GffCompare always loads the whole transcript data in memory for clustering and other analysis). One may only be interested to know _if_ and _how_ these many transcripts overlap the reference annotation, and further analyze only those which have specific types of overlap with the reference annotation transcripts (or none at all, i.e. if they do not overlap any of it, which may be the case for putative _novel_ transcripts). 
That's where the `trmap` utility comes in, as this program reports, for each query transcript, all the reference overlaps found, along with their _classification codes_ as described in the GffCompare documentation. The main feature of 'trmap' is that it allows _streaming_ of a very large file of query transcripts (in GFF or BED format) to be checked and classified against a reference annotation file (again, in GFF or BED format).

`trmap` ("Tanscript vs. Refererence MAPping") first loads the reference annotation file in memory as an interval tree and then streams the query file (which can be also provided at stdin) while checking and reporting any overlaps found, and classifies the relationship with reference transcripts using a (subset of) the "class codes" like those assigned by gffcompare (see http://ccb.jhu.edu/software/stringtie/gffcompare.shtml). 

The streaming input GFF query input file to be streamed must be _well-formed_ -- i.e. exons MUST be grouped together by transcript ID and immediately follow their parent feature if present. (for BED this is always the case due to the fact that exons are embedded in the same line).

## Building from source
Steps for building this package from github:
```
  cd /some/build/dir
  git clone https://github.com/gpertea/gffcompare
  cd gffcompare
  make release
```
If you downloaded the standalone source package `gffread-*.tar.gz` then just unpack that, change to the unpacked directory and run `make release` there.


This should build the **gffcompare** and **trmap** binary executables in the 
current directory.

