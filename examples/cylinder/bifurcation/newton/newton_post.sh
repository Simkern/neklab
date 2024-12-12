ifile=./lightkrylov_nwt.log
grep -E "rnorm" $ifile > rnorm.txt
grep -E "linear solver converged" $ifile > gmres_step.txt
grep -E "converged|outer" $ifile | awk '{ print $10, " ", $12}' | sed 's/[a-z]//g' > gmres_step_outer.txt