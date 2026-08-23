#!/bin/bash
set -euo pipefail

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Directory for the release
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"
RELEASE_DIR="$SRC_DIR/../../releases"
VERSION_DIR="$RELEASE_DIR/singer-$VERSION-beta-linux-x86_64"
BUILD_DIR="$SRC_DIR/build"
JOBS=$(( $(nproc) > 32 ? 32 : $(nproc) ))

# Create version directory
rm -rf "$VERSION_DIR" "$BUILD_DIR"
mkdir -p "$VERSION_DIR" "$BUILD_DIR/release" "$BUILD_DIR/debug"

cd "$SRC_DIR"

SRC_DIRS="ARG HMM moves sampler utils"
INCLUDES="-I."
for d in $SRC_DIRS; do INCLUDES="$INCLUDES -I$d"; done
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

for d in $SRC_DIRS; do mkdir -p "$BUILD_DIR/release/$d" "$BUILD_DIR/debug/$d"; done

printf '%s\n' $SOURCES | xargs -P "$JOBS" -I CPPFILE \
    g++ -std=c++17 -O3 -g -flto=8 $INCLUDES -c CPPFILE -o "$BUILD_DIR/release/CPPFILE.o"
g++ -std=c++17 -O3 -g -flto=8 $(find "$BUILD_DIR/release" -name '*.o') -o "$VERSION_DIR/singer"

printf '%s\n' $SOURCES | xargs -P "$JOBS" -I CPPFILE \
    g++ -std=c++17 -g $INCLUDES -c CPPFILE -o "$BUILD_DIR/debug/CPPFILE.o"
g++ -std=c++17 -g $(find "$BUILD_DIR/debug" -name '*.o') -o "$VERSION_DIR/singer_debug"

# Copy additional files
cp singer_master "$VERSION_DIR/singer_master"
cp convert_to_tskit.py "$VERSION_DIR/convert_to_tskit.py"
cp index_vcf.py "$VERSION_DIR/index_vcf.py"
cp merge_ARG.py "$VERSION_DIR/merge_ARG.py"
cp convert_long_ARG.py "$VERSION_DIR/convert_long_ARG.py"
cp "$SRC_DIR/../../LICENSE" "$VERSION_DIR/LICENSE"
cp "$VERSION_DIR/singer" singer.tmp && mv -f singer.tmp singer
cp "$VERSION_DIR/singer_debug" singer_debug.tmp && mv -f singer_debug.tmp singer_debug

# Change directory to releases
cd "$RELEASE_DIR"

# Create a tarball with the version number in the name
rm -f "singer-$VERSION-beta-linux-x86_64.tar.gz"
tar -czf "singer-$VERSION-beta-linux-x86_64.tar.gz" "singer-$VERSION-beta-linux-x86_64"
echo "wrote $(pwd)/singer-$VERSION-beta-linux-x86_64.tar.gz ($(du -h "singer-$VERSION-beta-linux-x86_64.tar.gz" | cut -f1))"
#rm -rf "singer-$VERSION-beta-linux-x86_64"
