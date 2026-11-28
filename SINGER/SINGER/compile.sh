#!/bin/bash
set -euo pipefail

if [ "$#" -gt 1 ]; then
    echo "Usage: bash $0 [version_number]"
    exit 1
fi
VERSION=${1:-}

SRC_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SRC_DIR"

case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) PLATFORM=linux-x86_64; CXX=g++ ;;
    Darwin-arm64) PLATFORM=mac-arm64; CXX=clang++ ;;
    *) echo "Error: no build settings for $(uname -s) $(uname -m)."; exit 1 ;;
esac

SRC_DIRS="ARG HMM moves sampler utils"
INCLUDES="-I."
for d in $SRC_DIRS; do INCLUDES="$INCLUDES -I$d"; done
INCLUDES="$INCLUDES $(pkg-config --cflags htslib)"
if [ "$PLATFORM" = linux-x86_64 ]; then
    HTS_LIBS=$(pkg-config --static --libs-only-l htslib | tr ' ' '\n' | grep -v -x -e -lm -e -lpthread -e -pthread | tr '\n' ' ')
    TCMALLOC=""
    if pkg-config --exists libtcmalloc_minimal; then TCMALLOC="$(pkg-config --libs-only-L libtcmalloc_minimal) -ltcmalloc_minimal"; fi
    LIBS="$(pkg-config --libs-only-L htslib) -Wl,-Bstatic $HTS_LIBS -lz $TCMALLOC -Wl,-Bdynamic -lpthread -lm"
    RELEASE_FLAGS="-O3 -g -DNDEBUG -flto=8 -fno-math-errno"
    JOBS=$(( $(nproc) > 32 ? 32 : $(nproc) ))
else
    LIBS="$(pkg-config --libs htslib) -lz"
    RELEASE_FLAGS="-O3 -g -DNDEBUG -flto"
    JOBS=$(sysctl -n hw.ncpu)
fi

BUILD_DIR="$SRC_DIR/build"
rm -rf "$BUILD_DIR"
SOURCES=$(find . -path ./lab -prune -o -name '*.cpp' -print)

build() {
    local name=$1 flags=$2
    for d in $SRC_DIRS; do mkdir -p "$BUILD_DIR/$name/$d"; done
    printf '%s\n' $SOURCES | xargs -P "$JOBS" -I CPPFILE \
        $CXX -std=c++17 $flags $INCLUDES -c CPPFILE -o "$BUILD_DIR/$name/CPPFILE.o"
    $CXX -std=c++17 $flags $(find "$BUILD_DIR/$name" -name '*.o') $LIBS -o "$BUILD_DIR/$name/$name"
    mv -f "$BUILD_DIR/$name/$name" "$name"
}

build singer "$RELEASE_FLAGS"
build singer_debug "-g"
rm -rf "$BUILD_DIR"
echo "wrote $SRC_DIR/singer and $SRC_DIR/singer_debug"

if [ -n "$VERSION" ]; then
    NAME="singer-$VERSION-beta-$PLATFORM"
    RELEASE_DIR="$SRC_DIR/../../releases"
    rm -rf "$RELEASE_DIR/$NAME"
    mkdir -p "$RELEASE_DIR/$NAME"
    cp singer singer_debug singer_master python_scripts/*.py python_scripts/*.pyx "$SRC_DIR/../../LICENSE" "$RELEASE_DIR/$NAME/"
    cd "$RELEASE_DIR"
    rm -f "$NAME.tar.gz"
    tar -czf "$NAME.tar.gz" "$NAME"
    echo "wrote $(pwd)/$NAME.tar.gz ($(du -h "$NAME.tar.gz" | cut -f1))"
    rm -rf "$NAME"
fi
