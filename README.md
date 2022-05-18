# Pion effiency analysis

Provides a `DSelector` and tools for running it for studying the charged pion
efficiency at the GlueX experiment.

## How to run

The `DSelector` can be run on the JLab cluster, or locally. Whatever method is
used, make sure to first source the environment file in the `env` directory.

### Locally

Use ROOT to run the `dselector/run_simple.C` with a config file from the
`dselector/config` directory as a parameter. The config file should consist of
an analysis tree name and file name, separated by spaces.

### Cluster

Use the `job/launch-*.csh` scripts to submit jobs to the cluster. To produce
analysis trees using `hd_root`, use the `job/launch-hdroot-*.csh` scripts. To
run the `DSelector` on the analysis trees, use the `job/launch-dselector-*.csh`
scripts. The scripts take no arguments. For running on different datasets, copy
and modify these launch scripts, including a configuration file (located in
`job/config-hdroot` or `job/config-dselector`) and run numbers.

The `job/env` folder contains scripts that set up the environment needed by each
of the jobs. Specifically, the environment variable `ACCID_METHOD` is used by
the `DSelector` to correctly apply accidental subtraction.

The `job/jana` folder contains various jana configuration files.

## Output

The `DSelector` produces a series of histograms organized into folders by the
types of cuts that were applied. Some naming conventions used are:

* `Pi0` refers to the neutral pion (both photons from the decay are measured).
* `Pi1` refers to the "found" charged pion, which is always detected.
* `Pi2` refers to the "missing" charged pion, which is only sometimes detected.
* `Found` suffix for histograms where both pions are found.
* `Missing` suffix for histograms where only one pion was found.
* `Measured` suffix for histograms using measured momenta only.
* `Mixed` suffix for histograms using mixture of measured and kinematic-fitted
  momenta.
* `MMOP` refers to the "Missing Mass Off Proton" method of estimating omega
  yield.
* `M3PI` refers to the "Mass of 3 Pions" method of estimating omega yield.

Unless otherwise specified in the name, histograms contain both `Found` and
`Missing` events, and use only the kinematic-fitted momenta.

## Modifying the code

The code can be found in the `dselector/DSelector_omega_pi_effic.C` file and its
corresponding header. It is mostly similar to the default `DSelector` structure.
However, there are a few aspects which need some additional explanation.

### Cuts

Cuts are defined by the `CUT_VARIANT` variable defined in the header file. Each
row in this array corresponds to a different combination of cuts which will be
used when producing histograms. The rows have three parts:
* A name describing the cut.
* A series of flags for the "standard" cuts, which are applied to all events.
* A series of flags for the "found" cuts, which are applied only to events in
  which both pions are found.

The `CUT_VARIANT` variable only determines which cuts are included in any given
histogram. The cuts themselves are defined in the `locCut[]` and `locCutFound[]`
variables from the source file. Each flag in the `CUT_VARIANT` variable refers
to whether the cut at a specific index in the `locCut[]` and `locCutFound[]`
variables should be applied.

To add a new combination of cuts, simply add a new row to the `CUT_VARIANT`
variable.

To add a new type of cut, first increase the `CUT_COUNT` or `CUT_FOUND_COUNT`
variable in the header. Then, add the new flag for the cut to the `CUT_VARIANT`
variable. Finally, add the calculation for the cut to the `locCut[]` or
`locCutFound[]` variables.

### Uniqueness tracking

Uniqueness tracking is done by the `UniquenessTracker` type defined in the
source file. Each time a combo needs to be checked for uniqueness, the `Check`
method is used, which automatically updates the internal map used for uniqueness
tracking. This type is equivalent to the normal way of doing uniqueness
tracking.

### Histograms

When adding a new histogram, it should be declared as an array of size
`CUT_VARIANT_COUNT`, as done with the existing histograms in the header. The
histogram must be defined in the loop labeled `CREATE HISTOGRAMS`. The
histograms don't need unique names for each cut combination, because all
histograms for a specific cut combination are placed in their own file. Then,
the histogram should be filled in one of the two loops labeled as
`FILL HISTOGRAMS`. One loop is for "found" events only, and the other is for all
events. These loops automatically apply all of the cuts corresponding to any cut
combination, so that only the surviving events will be included in the new
histogram.

