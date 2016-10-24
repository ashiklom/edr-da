#!/bin/bash

dbh=$1
pft=$2
dens=${3:-0.05}

pushd run-ed/1cohort/dens${dens}/dbh${dbh}/${pft}
./ed_2.1
popd
