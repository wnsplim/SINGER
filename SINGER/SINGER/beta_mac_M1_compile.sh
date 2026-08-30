#!/bin/bash
set -euo pipefail

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1
VERSION_DIR="../../releases/singer-$VERSION-beta-mac-arm64"
rm -rf "$VERSION_DIR"
mkdir -p "$VERSION_DIR"

# Compile the program with optimizations and debugging information
INCLUDES="-I. -IARG -IHMM -Imoves -Isampler -Iutils"
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

clang++ -std=c++17 -O3 -g -DNDEBUG -flto $INCLUDES $SOURCES -o "$VERSION_DIR/singer"

# Compile the debug version of the program
clang++ -std=c++17 -g $INCLUDES $SOURCES -o "$VERSION_DIR/singer_debug"

# Copy additional files
cp singer_master "$VERSION_DIR/singer_master"
cp convert_to_tskit.py "$VERSION_DIR/convert_to_tskit.py"
cp index_vcf.py "$VERSION_DIR/index_vcf.py"
cp merge_ARG.py "$VERSION_DIR/merge_ARG.py"
cp convert_long_ARG.py "$VERSION_DIR/convert_long_ARG.py"
cp ../../LICENSE "$VERSION_DIR/LICENSE"
cp "$VERSION_DIR/singer" singer.tmp && mv -f singer.tmp singer
cp "$VERSION_DIR/singer_debug" singer_debug.tmp && mv -f singer_debug.tmp singer_debug

# Change directory to releases
cd ../../releases

# Create a tarball with the version number in the name
# Rename to indicate ARM architecture (macOS M1)
rm -f "singer-$VERSION-beta-mac-arm64.tar.gz"
tar -czf "singer-$VERSION-beta-mac-arm64.tar.gz" "singer-$VERSION-beta-mac-arm64"
echo "wrote $(pwd)/singer-$VERSION-beta-mac-arm64.tar.gz ($(du -h "singer-$VERSION-beta-mac-arm64.tar.gz" | cut -f1))"
rm -rf "singer-$VERSION-beta-mac-arm64"
