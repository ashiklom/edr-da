#!/bin/bash

cohorts=$1
dbh=$2
pft=$3
dens=${4:-0.05}

pushd run-ed/${cohorts}/dens${dens}/dbh${dbh}/${pft}
#./ed_2.1
./ed_2.1-opt
popd
