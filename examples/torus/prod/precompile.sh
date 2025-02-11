#!bash
source ~/.bashrc

cd compile
for cdir in *_*; do
   cd $cdir
   if [ ! -f "torus.usr" ]; then
      cp -v torus.usr_ref torus.usr
   fi
   makeneklab
   cd -
done
cd ..

pwd              
