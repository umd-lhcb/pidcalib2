{
  description = "pidcalib2.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
  };

  outputs = { self, root-curated, nixpkgs, flake-utils }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ root-curated.overlay self.overlay ];
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
            root
            # Python stack
            pidcalib2
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              gyre-fonts
            ];
          };

          shellHook = ''
            export MPLBACKEND=agg  # the backend w/o a UI
            export MPLCONFIGDIR=$(pwd)/.matplotlib

            # from https://discourse.nixos.org/t/python-package-with-runtime-dependencies/5522/2
            # still mess up with system tools with the libc are diifferent
            #export LD_LIBRARY_PATH=${pkgs.lib.makeLibraryPath [pkgs.stdenv.cc.cc]}
          '';
        };
      });
}
