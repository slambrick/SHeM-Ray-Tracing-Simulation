#!/bin/bash

make bin/simple
make bin/random

echo Limit of the gaussian tail?
read limit

./bin/simple
./bin/random $limit

./plot_dist.R $limit
