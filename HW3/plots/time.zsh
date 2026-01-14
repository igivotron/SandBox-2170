#!/bin/zsh

for x in {50..400..10}; do
    python main.py -N $x
done