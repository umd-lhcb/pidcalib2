{
  description = "pidcalib2.";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
        };
        devShell = pkgs.mkShell rec {
          name = "pidcalib2-dev";
          buildInputs = (with pkgs; with pythonPackages; [
            # Python stack
            virtualenvwrapper
            boost-histogram
            logzero
            matplotlib
            numpy
            pandas
            tqdm
          ]);

          shellHook = ''
            export PATH=$(pwd)/scripts:$(pwd)/bin:$PATH

            # Allow the use of wheels.
            SOURCE_DATE_EPOCH=$(date +%s)

            VENV=./.virtualenv

            if test ! -d $VENV; then
              virtualenv $VENV
            fi
            source $VENV/bin/activate

            # allow for the environment to pick up packages installed with virtualenv
            export PYTHONPATH=$VENV/${python.sitePackages}/:$PYTHONPATH
            # from https://discourse.nixos.org/t/python-package-with-runtime-dependencies/5522/2
            # still mess up with system tools with the libc are diifferent
            #export LD_LIBRARY_PATH=${pkgs.lib.makeLibraryPath [pkgs.stdenv.cc.cc]}
          '';
        };
      });
}
