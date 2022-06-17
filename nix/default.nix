{ stdenv
, buildPythonPackage
, boost-histogram
, logzero
, matplotlib
, mplhep
, numpy
, pandas
, tqdm
, uproot
}:

# FIXME: We require uproot 4 but it's not in nixpkgs yet.

buildPythonPackage rec {
  pname = "pidcalib2";
  version = "1.0.6";

  src = builtins.path { path = ./..; name = pname; };

  propagatedBuildInputs = [
    boost-histogram
    logzero
    matplotlib
    mplhep
    numpy
    pandas
    tqdm
    uproot
  ];

  doCheck = false;
}
