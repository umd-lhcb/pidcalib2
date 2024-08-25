# pidcalib2

This README contains instructions specific for UMD LHCb group.

For the project's upstream README, see [here](https://gitlab.cern.ch/lhcb-rta/pidcalib2).


## Usage

1. Make sure you are on `glacier`, our server
2. Clone this project
3. Type `nix develop` in the project root (if you don't have `direnv` set up).
4. Take a look at `efficiency_gen` folder, where some example scripts are stored

## Resave PIDCalib ntuples so that they can be opened in ROOT 5

We are only able to [calculate the `uBDT`](https://github.com/umd-lhcb/MuonBDTPid) in ROOT 5. The latest
PIDCalib ntuples for 2017 and 2018 were generated in ROOT 6.30.04, and are not backwards-compatible
with ROOT 5, but they are openable in ROOT 6.24. Thus, to use them in ROOT 5 they need to be resaved
in 6.24 with the following script ini `glacier` where we copied them to

```
nix develop
make
./scripts/resave_all_pidcalib_ntuples.py -i /home/public/pidcalib_ntuples/remote -t 17 -o /home/public/pidcalib_ntuples/remote/resaved/
./scripts/resave_all_pidcalib_ntuples.py -i /home/public/pidcalib_ntuples/remote -t 18 -o /home/public/pidcalib_ntuples/remote/resaved/
```
