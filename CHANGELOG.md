ChangeLog
=========

A simple change log file to keep track on what is new.


Versions
========

* 0.6
    * new commandline parameter: `--min-substantive-abundance` (a very important new parameter for better noise control, see http://oligotyping.org/2012/09/18/command-line-parameters-explained/ for details).
    * new commandline parameter: `--gen-oligotype-sets` (oligotype sets are going to be generated only when requested).
    * new commandline parameter: `--exclude-oligotypes`.
    * naming convention for agglomerated oligotypes has been changed; now they are known as oligotype sets.

* 0.5
    * grouping oligotypes based on their frequency patterns across datasets has been added to the pipeline.
    * code directory structure has been changed.

* 0.4
	* visualization showing normalized oligotype distribution across datasets has been added to the pipeline.

* 0.3
	* new commandline parameter: `--colors-list-file`.

* 0.2
	* unittests implemented for various functions
	* class structure was improved

* 0.1
	* initial oligotyping library
