{ nixpkgs ? (import (builtins.fetchTarball {
  name = "matthuszagh-nixpkgs-pyspice-2020-05-23";
  url = "https://github.com/matthuszagh/nixpkgs/archive/8516c750b95b0b97538445df84f264ba67792834.tar.gz";
  sha256 = "0fhg3878dj0ql7nzsmqawnf1q0xr22zabch977ffvi647ly3ga64";
}) {})
}:

let
  custompkgs = import <custompkgs> {};
  pkgs = nixpkgs;
  pythonEnv = (pkgs.python3Full.buildEnv.override {
    extraLibs = with pkgs.python3Packages; [
      numpy
      matplotlib
      pathos
      scipy
    ] ++ (with custompkgs; [
      circlib
    ]);
    ignoreCollisions = true;
  });
in
pkgs.mkShell rec {
  buildInputs = with pkgs; [
    pythonEnv
    qucs
  ];

  KICAD_SYMBOL_DIR="/home/matt/src/kicad-symbols";
}
