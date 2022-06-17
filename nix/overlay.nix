final: prev:

{
  pythonOverrides = prev.lib.composeExtensions prev.pythonOverrides (finalPy: prevPy: {
    pidcalib2 = finalPy.callPackage ./default.nix { };
    tqdm = finalPy.callPackage ./tqdm { };
  });
  python3 = prev.python3.override { packageOverrides = final.pythonOverrides; };
}
