#!/bin/bash

INCLUDES="-I. -IARG -IHMM -Imoves -Isampler -Iutils"
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

g++ -std=c++17 -O3 -g -DNDEBUG -static -flto=8 -fno-math-errno $INCLUDES $SOURCES -o singer
g++ -std=c++17 -g -static $INCLUDES $SOURCES -o singer_debug

