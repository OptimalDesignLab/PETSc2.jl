# download petsc ahead of time if needed

petsc_name="petsc-3.6.0"
fmt=".tar.gz"
url="http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/"

url_complete = string(url, petsc_name, fmt)
download(url_complete)


