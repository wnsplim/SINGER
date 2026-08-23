#!/bin/bash

INCLUDES="-I. -IARG -IHMM -Imoves -Isampler -Iutils"
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

g++ -std=c++17 -O3 -g -static $INCLUDES $SOURCES -o singer
g++ -std=c++17 -g -static $INCLUDES $SOURCES -o singer_debug

