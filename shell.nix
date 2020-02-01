{ nixpkgs ? (import <nixpkgs> {})
}:

let
  custompkgs = import <custompkgs> {};
  pkgs = (nixpkgs // custompkgs);
  libcircuit = pkgs.libcircuit;
  mh-python = pkgs.python3.withPackages (ps: with ps; [
    numpy
    matplotlib
    pathos
    scipy
  ] ++ (with custompkgs; [
    libcircuit
    skidl
  ]));
  kicad = pkgs.kicad;
in
pkgs.mkShell rec {
  buildInputs = with pkgs; [
    python3Full
    mh-python
    python-openems
    python-csxcad

    # pcb cad
    kicad

    # ems
    (openems.override {withMPI = false; })
    appcsxcad
    hyp2mat
  ];

  KICAD_SYMBOL_DIR="${kicad.out}/share/kicad/library";
}
