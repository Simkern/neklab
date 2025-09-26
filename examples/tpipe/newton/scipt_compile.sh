RUN=`pwd`
cp tpipe.usr compile
cd compile
../../../app/makeneklab
cd $RUN
cp -v compile/nek5000 .
