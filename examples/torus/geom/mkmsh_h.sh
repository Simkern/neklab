if [ $# -ne 1 ]; then echo -e "\nWrong argument. Abort.\n"; exit 1; fi
#
#  Usage: bash mkmsh_h.sh torus_coarse2Dh
#
CNAME=$1
TDIR="templates"
LOGFILE="mkms_h.log"
if [ -e $LOGFILE ]; then rm $LOGFILE; fi

echo "$CNAME:"
echo "  n2to3 ..."
## n2to3
sed "s/casename2D/${CNAME}/g" ${TDIR}/n2to3_cmd_template.txt > n2to3_cmd.txt;
sed -i "s/casename3D/${CNAME}_/g" n2to3_cmd.txt;
sed -i "s/2D_/3D/g" n2to3_cmd.txt;
n2to3 < n2to3_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  genmap ..."
## genmap
echo "${CNAME}_" > genmap_cmd.txt
echo "0.000001" >> genmap_cmd.txt
sed -i "s/2D_$/3D/g" genmap_cmd.txt
genmap < genmap_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "done."
# cleanup
rm -f fort.99 gmsh2nek_cmd.txt n2to3_cmd.txt re2torea_cmd.txt genmap_cmd.txt
