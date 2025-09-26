# T-pipe example

Files:
├── compile
│   └── SIZE                                    ! definitions of static arrays dimensions
├── geom
│   └── tpipe.re2                               ! mesh file binary
├── newton
│   ├── compile -> ../compile
│   ├── tpipe.re2 -> ../geom/tpipe.re2
│   ├── scipt_compile.sh                        ! compile script
│   ├── tpipe.par                               ! Nek5000 parameter file
│   └── tpipe.usr                               ! Nek5000 setup source file
└── README.md

To compile the code:
* ensure that `Nek5000` has been successfully cloned and the `LightKrylov`-specific changes have been executed. This is most easily acheived by running the `Nek5000_setup.sh` script in the neklab root directory (Note: the script must be executable).
* ensure that `LightKrylov` has been successfully cloned and installed on the system. This is most easily acheived by running the `LightKrylov_setup.sh` script in the neklab root directory (Note: the script must be executable). We recommend running the test suite to check that everything works correctly.
* build the case using script `script_compile`.

To run the case:
* generate the processor map file `tpipe.ma2` using the tool `genmap` distributed together with `Nek5000`.
* run the code using the parallel executable provided in `Nek5000/bin/`.