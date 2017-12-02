#!/bin/bash

# This script must be run from the /docs directory
# It will build the documentation using Documenter and then deploy it to
# github 
# $1 is the version name of the docs (for example, latest, v0.1, v0.2 etc.)
# defaults to "latest" if not specified
# figure out the name of the branch/tag to call the docs

echo "number of args = $#"
if [[ $# > 1 ]]; then
  echo "Usage: deploy_docs.sh doc_name"
  exit 1
fi

if [[ $# == 0 ]]; then
  dname="latest"
else
  dname=$1
fi

echo "deploying docs named: $dname"

start_dir=`pwd`
julia ./make.jl

if [ $? != 0 ]; then
  echo "Error: Building docs failed, aborting ..."
  exit 1
fi

commit_msg=`git log -n 1 --abbrev-commit --oneline`

# this incantation should work on OS X and Linux
tdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tdir'`
echo "using temporary directory = $tdir"

# download another copy of the repo to a temporary directory
cd $tdir

if [ $? != 0 ]; then
  echo "Error: could not cd into temporary directory"
  exit 1
fi

git clone https://github.com/OptimalDesignLab/PETSc2.jl.git
cd ./PETSc2.jl
git checkout gh-pages
mkdir -vp $dname  # make the directory if it doesn't exist
cd ./$dname
rm -r ./*  # remove everything

# copy docs to temporary repo
pres_dir=`pwd`
echo "copying files from $start_dir/build to $pres_dir"
cp -r $start_dir/build/* $pres_dir
git add ./*
git commit -am "building docs from: $commit_msg"

# deploy
echo "deploying to github"
git push origin gh-pages --force

cd $start_dir

rm -rf $tdir  # clean up temporary directory
