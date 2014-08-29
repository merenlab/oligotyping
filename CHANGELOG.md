ChangeLog
=========

A simple change log file to keep track on what is new.


Versions
========

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
