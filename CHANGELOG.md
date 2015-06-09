ChangeLog
=========

A simple change log file to keep track on what is new.


Versions
========

* 2.1 (NOT RELASED YET)
    * Proper handling of files with 0 entropy.

* 2.0 (2015-05-22)
    * Fixing a ancient bug: remove carriage return characters from sample mapping files (so people who export their mapping files from EXCEL on Mac computers can get nice figures as well).

* 1.9 (2015-03-31)
    * New script `o-smart-trim` to trim ragged ends of alignments with care.
    * Slightly better command line interface.
    * Minor bug fixes.

* 1.7 (2015-01-20)
    * Incompatibility with pip 6.0.6 is fixed.
    * Minor bug fixes.

* 1.6 (2014-11-30)
    * Minor bug fixes.

* 1.5 (2014-11-11)
    * o-dist.R is removed from the codebase. Another script that computes pairwise distances of sequencs in a given FASTA file to generate a distance matrix is added (o-dist) along with another script to visualize this informatopn (o-visualize-distance-matrix.R).

* 1.4 (2014-10-01)
    * A minor update to remove a once-useful command line parameter from the oligotyping pipeline, --gen-sample-oligo-networks. The same functionality is now available as a separate script (o-generate-oligo-base-networks). I thank Sinisa Bratulic for reporting it.

* 1.3
    * A minor update to fix the broken backward compatibility of Django v1.7. I thank Bernd from Amsterdam for pointing it out.

* 1.2
    * Bug in MANIFEST file is fixed.
    * Bug related to --quick flag (which resulted in under-estimaded benchmark scores) is fixed.

* 1.1
    * ggplot2-related bug (which resulted in blank stackbar plots) is fixed.

* 1.0
    * Minimum Entropy Decomposition is an part of the pipeline.
    * Experimental purity score implemented by Doğancan Özturan is now included in the HTML output.

* 0.96
    * Improved HTML output.
    * Better error reporting.

* 0.95
    * GEXF support. Now oligotyping pipeline generates Gephi (https://gephi.org/) compatible XML file for network analysis.
    * Bug fixes. 

* 0.9
    * Multithreading support. Oligotyping pipeline runs processes in parallel whenever it is possible.
    * More comprehensible HTML output.
    * RUNINFO.cPickle format has changed.

* 0.8
    * Dependency to USEARCH is removed from the pipeline. NCBI's blastn is now the default search tool for relevant operations.
    * Reads that are being discarded during the noise removal step are now tracked to avoid biases.
    * New commandline parapmeter: `--sample-mapping` (which supports TAB delimited categorical mapping files for samples to be used for post-analysis visualizations; a relevant blog post: http://oligotyping.org/2013/02/04/basic-sample-mapping-with-oligotyping/).

* 0.7
    * Major performance improvements.
    * All output matrices are transposed.

* 0.6
    * New commandline parameter: `--min-substantive-abundance` (a very important new parameter for better noise control, see http://oligotyping.org/2012/09/18/command-line-parameters-explained/ for details).
    * New commandline parameter: `--gen-oligotype-sets` (oligotype sets are going to be generated only when requested).
    * New commandline parameter: `--exclude-oligotypes`.
    * Naming convention for agglomerated oligotypes has been changed; now they are known as oligotype sets.

* 0.5
    * Grouping oligotypes based on their frequency patterns across datasets has been added to the pipeline.
    * Code directory structure has been changed.

* 0.4
	* Visualization showing normalized oligotype distribution across datasets has been added to the pipeline.

* 0.3
	* New commandline parameter: `--colors-list-file`.

* 0.2
	* Unittests implemented for various functions
	* Class structure was improved

* 0.1
	* Initial oligotyping library
