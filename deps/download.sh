#!/bin/bash
# download petsc ahead of time if needed

petsc_name="petsc-3.9.4"
fmt=".tar.gz"
url="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/"

wget "$url""$petsc_name""$fmt"


