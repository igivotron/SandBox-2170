#!/bin/zsh

for x in {50..400..10}; do
    for y in {1..10..1}; do
        python ./poisson.py -M $x
    done
done