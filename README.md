Oligotyping
===============================================================================

Oligotyping offers a means to utilize extremely subtle variants in amplicon
data using information theory without falling into the trap of pairwise
sequence alignments and counting mismatches to estimate distances between
reads. This approach helped a lot to recove very subtle nucleotide variation
among 16S Ribosomal RNA gene amplicons when the only available alternative for
amplicon analyses were 3% OTUs or taxon assignments. Oligotyping and MED revealed
subtle but critical ecological patterns, and hidden diversity concealed within
environmental datasets.

The application of entropy-based decomposition of data is not limited to amplicon
data and can also be used in metagenomic contexts. But admittedly with the
shifting priorities of Meren Lab from amplicon-based analyses to metagenomics
slowed down the development of the oligotpying pipeline.

If you are still planning to use this resource, get in touch with Meren so he can
help you get started with some useful best practices material.

Installation is still possible, and it -most likely- still works well :)

Installation
===============================================================================

The recommended way to install oligotyping is within an isolated
[conda](https://docs.conda.io/en/latest/miniconda.html) environment. This keeps
the Python packages, the BLAST tools, and the R libraries the pipeline relies on
together in one place, and away from the rest of your system.

Oligotyping runs on Python 3.7 and later, and is tested against current stable
releases (Python 3.13 / 3.14), so it is safe to create the environment with the
latest Python your conda channels provide:

```bash
# create and activate a fresh environment (any Python >= 3.7 works;
# `python=3` picks the latest stable release available to you)
conda create -y -n oligotyping python=3
conda activate oligotyping

# non-Python dependencies:
#   - BLAST (blastn / makeblastdb) is required for both the oligotyping and the
#     minimum entropy decomposition (MED) pipelines
#   - R (with the packages below) is used to generate the figures
conda install -y -c bioconda blast
conda install -y -c conda-forge \
    r-base r-vegan r-ggplot2 r-optparse r-pheatmap \
    r-mass r-gtools r-reshape r-rcolorbrewer r-compute.es

# get the latest development version from source and install it in editable mode
mkdir -p ~/github && cd ~/github
git clone https://github.com/merenlab/oligotyping.git
cd oligotyping
pip install -e .
```

Once the installation is complete, confirm that the entry points are available
(remember to `conda activate oligotyping` in every new terminal session first):

```bash
oligotype --version
decompose --version
entropy-analysis --version
```


Notes for macOS users
-------------------------------------------------------------------------------

A couple of macOS-specific things are worth keeping in mind:

* **Xcode Command Line Tools.** Cloning the repository (and building any dependency
  that is not distributed as a pre-built wheel) requires Apple's Command Line Tools.
  If `git` is missing, or you run into compiler errors during installation, install
  them with `xcode-select --install`. If you have the full Xcode installed but have
  never accepted its license, you may additionally need to run `sudo xcodebuild -license`.

* **matplotlib backend.** oligotyping renders some of its figures from within worker
  processes while running multi-threaded, and the interactive backend that matplotlib
  tends to select by default on macOS (`macosx`) is not safe to use that way. If you
  run into figure-related crashes or hangs, force the non-interactive `Agg` backend by
  exporting an environment variable before you run the pipeline:

  ```bash
  export MPLBACKEND=Agg
  ```

  You can add that line to your shell profile (or to the activation hook below) so it
  is always set within the environment.


Following the latest development version
-------------------------------------------------------------------------------

Because the step above installs oligotyping in "editable" mode (`pip install -e .`),
the code that actually runs is always whatever is in your local clone at
`~/github/oligotyping`. You can make that clone follow the latest development code
automatically by adding a small activation hook, so that the repository is synchronized
with GitHub every time you activate the environment:

```bash
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d/
cat <<EOF >${CONDA_PREFIX}/etc/conda/activate.d/oligotyping.sh
# every time this conda environment is activated, synchronize the local
# oligotyping clone with the latest code from the GitHub repository:
echo -e "\033[1;34mUpdating oligotyping from GitHub \033[0;31m(press CTRL+C to cancel)\033[0m ..."
cd ~/github/oligotyping && git pull && cd - > /dev/null
EOF
```

(Adjust the `~/github/oligotyping` path above if you cloned the repository elsewhere.)
From now on, each `conda activate oligotyping` will pull the newest changes before you
run anything.


Licence
===============================================================================

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
