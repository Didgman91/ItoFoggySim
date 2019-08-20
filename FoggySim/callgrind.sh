#!/bin/bash

outputLog="callgrind.log"
outputCallgrind="callgrind.profile"

SECONDS=0

#valgrind --tool=callgrind --log-file="${outputLog}" --callgrind-out-file="${outputCallgrind}" ./Debug/FoggySim
valgrind --tool=callgrind --collect-systime=yes --log-file="${outputLog}" --callgrind-out-file="${outputCallgrind}" ./Debug/FoggySim

echo $SECONDS

kcachegrind "${outputCallgrind}"
