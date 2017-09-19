#!/bin/bash
# download petsc ahead of time if needed

petsc_name="petsc-3.7.6"
fmt=".tar.gz"
url="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/"

wget "$url""$petsc_name""$fmt"


