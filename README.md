# PIDCalib2

A set of software tools for estimating LHCb PID efficiencies.

The package includes several user-callable modules:
- `make_eff_hists` creates histograms that can be used to estimate the PID efficiency of a user's sample
- `ref_calib` calculates the LHCb PID efficiency of a user reference sample
- `merge_trees` merges two ROOT files with compatible `TTree`s
- `pklhisto2root` converts [Pickled](https://docs.python.org/3/library/pickle.html) [boost-histogram](https://github.com/scikit-hep/boost-histogram)s to ROOT histograms

## Setup

When working on a computer that where the LHCb software stack is available (LXPLUS, some university cluster, etc.), one can setup PIDCalib2 by running
```sh
lb-conda pidcalib bash
```
After this, the following commands will be available
```sh
pidcalib2.make_eff_hists
pidcalib2.ref_calib
pidcalib2.merge_trees
pidcalib2.pklhisto2root
```
To run `make_eff_hists`, you will need access to CERN EOS. You don't need to do anything special on LXPLUS. On other machines, you will usually need to obtain a Kerberos ticket by running
```sh
kinit [username]@CERN.CH
```

### Installing from PyPI

The PIDCalib2 package is available on [PyPI](https://pypi.org/project/pidcalib2/). It can be installed on any computer via `pip` simply by running (preferably in a virtual environment; see [`venv`](https://docs.python.org/3/library/venv.html))
```sh
pip install pidcalib2
```
Note that this will install the [`xrootd`](https://pypi.org/project/xrootd/) *Python bindings*. One also has to install XRootD itself for the bindings to work. See [this page](https://xrootd.slac.stanford.edu/index.html) for XRootD releases and instructions.

## `make_eff_hists`

This module creates histograms that can be used to estimate the PID efficiency of a user's sample.

Reading all the relevant calibration files can take a long time. When running a configuration for the first time, we recommend using the `--max-files 1` option. This will limit PIDCalib2 to reading just a single calibration file. Such a test will reveal any problems with, e.g., missing variables quickly. Keep in mind that you might get a warning about empty bins in the total histogram as you are reading a small subset of the calibration data. For the purposes of a quick test, this warning can be safely ignored.

### Options

To get a usage message listing all the options, their descriptions, and default values, type
```
pidcalib2.make_eff_hists --help
```

The calibration files to be processed are determined by the `sample`, `magnet`, and `particle` options. All the valid combinations can be listed by running
```sh
pidcalib2.make_eff_hists --list configs
```

Aliases for standard variables are defined to simplify the commands. We recommend users use only the aliases when specifying variables. When you use a name that isn't an alias, a warning message like the following will show up in the log. Use with caution.
```
'probe_PIDK' is not a known PID variable alias, using raw variable
```
All aliases can be listed by running
```sh
pidcalib2.make_eff_hists --list aliases
```

A file with alternative binnings can be specified using `--binning-file`. The file must contain valid JSON specifying bin edges. For example, a two-bin binning for particle `Pi`, variable `P` can be defined as
```json
{"Pi": {"P": [10000, 15000, 30000]}}
```
An arbitrary number of binnings can be defined in a single file.

Complex cut expressions can be created by chaining simpler expressions using `&`. One can also use standard mathematical symbols, like `*`, `/`, `+`, `-`, `(`, `)`. Whitespace does not matter.

### Examples:
- Create a single efficiency histogram for a single PID cut
  ```sh
  pidcalib2.make_eff_hists --sample Turbo18 --magnet up --particle Pi --pid-cut "DLLK > 4" --bin-var P --bin-var ETA --bin-var nSPDhits --output-dir pidcalib_output
  ```

- Create multiple histograms in one run (most of the time is spent reading
in data, so specifying multiple cuts is much faster than running
make_eff_hists sequentially)
  ```sh
  pidcalib2.make_eff_hists --sample Turbo16 --magnet up --particle Pi --pid-cut "DLLK > 0" --pid-cut "DLLK > 4" --pid-cut "DLLK > 6" --bin-var P --bin-var ETA --bin-var nSPDhits --output-dir pidcalib_output
  ```

- Create a single efficiency histogram for a complex PID cut
  ```sh
  pidcalib2.make_eff_hists --sample Turbo18 --magnet up --particle Pi --pid-cut "MC15TuneV1_ProbNNp*(1-MC15TuneV1_ProbNNpi)*(1-MC15TuneV1_ProbNNk) < 0.5 & DLLK < 3" --cut "isMuon==0" --bin-var P --bin-var ETA --bin-var nSPDhits --output-dir pidcalib_output
  ```

## `ref_calib`

This module uses the histograms created by `make_eff_hists` to assign efficiency to events in a reference sample supplied by the user. Adding of efficiency to the user-supplied file requires PyROOT and is optional.

The module works in two steps:

1. Calculate the efficiency and save it as a TTree in a separate file.
2. Optionally copy the efficiency TTree to the reference file and make it a friend of the user's TTree. The user must request the step by specifying `--merge` on the command line.

Be aware that `--merge` will modify your file. Use with caution.

### Options

The `sample` and `magnet` options are used solely to select the correct PID efficiency histograms. They should therefore mirror the options used when running `make_eff_hists`.

`bin-vars` must be a dictionary that relates the binning variables (or aliases) used to make the efficiency histograms with the variables in the reference sample. We assume that the reference sample branch names have the format `[ParticleName]_[VariableName]`. E.g., `D0_K_calcETA`, corresponds to a particle named `D0_K` and variable `calcETA`. If the user wants to estimate PID efficiency of their sample using 1D binning, where `calcETA` corresponds to the `ETA` binning variable alias of the calibration sample, they should specify `--bin-vars '{"ETA": "calcETA"}'`.

`ref-pars` must be a dictionary of particles from the reference sample to apply cuts to. The keys represent the particle branch name prefix (`D0_K` in the previous example), and the values passed are a list containing particle type and PID cut, e.g. `'{"D0_K" : ["K", "DLLK > 4"], "D0_Pi" : ["Pi", "DLLK < 4"]}'`.

The `--merge` option will copy the PID efficiency tree to your input file and make the PID efficiency tree a "Friend" of your input tree. Then you can treat your input tree as if it had the PID efficiency branches itself. E.g., `input_tree->Draw("PIDCalibEff")` should work. ROOT's "Friend" mechanism is an efficient way to add branches from one tree to another. Take a look [here](https://root.cern.ch/root/htmldoc/guides/users-guide/Trees.html#example-3-adding-friends-to-trees) if you would like to know more.

### Examples:
- Evaluate efficiency of a single PID cut and save it to `user_ntuple_PID_eff.root` without adding it to `user_ntuple.root`
  ```sh
  python -m pidcalib2.ref_calib --sample Turbo18 --magnet up --ref-file data/user_ntuple.root --output-dir pidcalib_output --bin-vars '{"P": "mom", "ETA": "Eta", "nSPDHits": "nSPDhits"}' --ref-pars '{"Bach": ["K", "DLLK > 4"]}'
  ```
- Evaluate efficiency of a single PID cut and add it to the reference file `user_ntuple.root`
  ```sh
  python -m pidcalib2.ref_calib --sample Turbo18 --magnet up --ref-file data/user_ntuple.root --output-dir pidcalib_output --bin-vars '{"P": "mom", "ETA": "Eta", "nSPDHits": "nSPDhits"}' --ref-pars '{"Bach": ["K", "DLLK > 4"]}' --merge
  ```
- Evaluate efficiency of multiple PID cuts and add them to the reference file
  ```sh
  python -m pidcalib2.ref_calib --sample Turbo18 --magnet up --ref-file data/user_ntuple.root --output-dir pidcalib_output --bin-vars '{"P": "P", "ETA": "ETA", "nSPDHits": "nSPDHits"}' --ref-pars '{"Bach": ["K", "DLLK > 4"], "SPi": ["Pi", "DLLK < 0"]}' --merge
  ```


## Development

1. Clone the repository from [GitLab](https://gitlab.cern.ch/lhcb-rta/pidcalib2)
2. (Optional) Set up a virtual environment
   ```sh
   python3 -m venv .venv
   source .venv/bin/activate
   ```
3. Install pinned dependencies
   ```sh
   pip install -r requirements-dev.txt
   ```
4. Install `xrootd` (possibly manually; see this [issue](https://github.com/xrootd/xrootd/issues/1397))
5. Run the tests
   ```sh
   pytest
   ```
6. Run the modules
   ```sh
   python3 -m src.pidcalib2.make_eff_hists -h
   ```

### Tips

Certain tests can be excluded like this
  ```sh
  pytest -m "not xrootd"
  ```
See available tags in the `src/pidcalib2/tests/test_*.py` files.
