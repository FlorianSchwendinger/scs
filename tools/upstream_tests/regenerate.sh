#!/bin/bash
# Regenerate inst/tinytest/problem_results/*.RDS from the C test suite bundled
# with SCS (src/scs/test/problems/*.h).
#
# Run from the package root:
#     bash tools/upstream_tests/regenerate.sh
#
# Steps: copy src/scs into tools/upstream_tests/build/, apply the src/scs_mods
# copies exactly as src/Makevars does, build libscsdir.a with the
# flags src/Makevars uses, compile dump_tests.c as an R shared object against
# it and R's BLAS/LAPACK, run upstream's test functions inside an R process
# (every scs()/scs_solve() call is recorded), and convert the records to RDS
# files.
set -u
ROOT=$(pwd)
if [ ! -f DESCRIPTION ] || [ ! -d src/scs ]; then
  echo "run this from the package root (DESCRIPTION and src/scs not found in $ROOT)"; exit 1
fi
HERE=$ROOT/tools/upstream_tests
BUILD=$HERE/build
TREE=$BUILD/src/scs
OUT=$ROOT/inst/tinytest/problem_results

echo "== 1. copy of the bundled sources with src/scs_mods applied (as src/Makevars does) =="
rm -rf "$BUILD"; mkdir -p "$BUILD/src"
cp -R src/scs "$TREE" && rm -rf "$TREE/.git"
cp -f src/scs_mods/scs.mk "$TREE/" && cp -f src/scs_mods/cones.c "$TREE/src/" || exit 1

echo "== 2. libscsdir.a with the package's flags =="
CC=$(R CMD config CC)
RCFLAGS=$(R CMD config CFLAGS)
PIC=$(R CMD config CPICFLAGS)
RINC=$(R RHOME)/include
PKGC="-I$TREE/include -I$RINC -DUSE_LAPACK -DR_LANG"
(cd "$TREE" && make purge > /dev/null 2>&1; \
 make OPT="" CC="$CC" CFLAGS="$RCFLAGS $PKGC $PIC" AR="$(R CMD config AR)" RANLIB="$(R CMD config RANLIB)" \
      out/libscsdir.a > "$BUILD/libscs.log" 2>&1) || { tail -20 "$BUILD/libscs.log"; exit 1; }

echo "== 3. dump_tests.so =="
cd "$BUILD" || exit 1
cp "$HERE/dump_tests.c" .
PKG_CPPFLAGS="-I$TREE/include -I$TREE/linsys -I$TREE/test -DUSE_LAPACK -DR_LANG -DCTRLC=1 -DGPU_TRANSPOSE_MAT=1 -Wno-unused-function" \
PKG_LIBS="$TREE/out/libscsdir.a $(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS) $(R CMD config FLIBS)" \
R CMD SHLIB dump_tests.c -o dump_tests.so > shlib.log 2>&1 || { cat shlib.log; exit 1; }

echo "== 4. upstream tests inside R =="
mkdir -p dumps
Rscript --vanilla "$HERE/run_dump.R" "$TREE" "$BUILD/dump_tests.so" "$BUILD/dumps" > run.log 2>&1; rc=$?
tail -1 run.log
if [ $rc -ne 0 ]; then grep "UPSTREAM ASSERTION FAILED" run.log; exit 1; fi

echo "== 5. RDS files -> $OUT =="
rm -f "$OUT"/*.RDS
Rscript --vanilla "$HERE/dumps_to_rds.R" "$BUILD/dumps" "$OUT" || exit 1
