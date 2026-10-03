#!/bin/bash

cd "$(dirname "$0")"
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH:+$PKG_CONFIG_PATH:}/space/s1/david9456/0.Refs/tools_install/lib/pkgconfig"
INCLUDES="-I. -IARG -IHMM -Imoves -Isampler -Iutils $(pkg-config --cflags htslib)"
TCMALLOC=""
if pkg-config --exists libtcmalloc_minimal; then TCMALLOC="$(pkg-config --static --libs libtcmalloc_minimal)"; fi
LIBS="$(pkg-config --static --libs htslib) $TCMALLOC -lz"
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

g++ -std=c++17 -O3 -g -DNDEBUG -static -flto=8 -fno-math-errno $INCLUDES $SOURCES $LIBS -o singer
g++ -std=c++17 -g -static $INCLUDES $SOURCES $LIBS -o singer_debug

