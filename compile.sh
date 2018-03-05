#!/bin/bash

mpic++ -std=c++11 src/genHap.cpp src/reader.cpp src/gene.cpp src/geneticOperation.cpp src/chromosome.cpp src/geneticAlgorithm.cpp -O3 -o GenHap