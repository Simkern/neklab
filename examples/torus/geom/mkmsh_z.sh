if [ $# -ne 2 ]; then echo -e "\nWrong argument. Abort.\n"; exit 1; fi
#
#  Usage: bash mkmsh_z.sh torus_coarse2D nz
#
CNAME=$1
NZ=$2
if [ ! -e ${CNAME}.re2 ]; then echo -e "\nMesh file ${CNAME}.re2 not found. Abort.\n"; exit 1; fi
TDIR="templates"
LOGFILE="mkmsh_z.log"
if [ -e $LOGFILE ]; then rm $LOGFILE; fi

echo "  n2to3 ..."
## n2to3
sed "s/casename2D/${CNAME}/g" ${TDIR}/n2to3_cmd_template.txt > n2to3_cmd.txt;
sed -i "s/casename3D/${CNAME}_${NZ}z_3D/g" n2to3_cmd.txt;
sed -i "s/2D_//g" n2to3_cmd.txt;
sed -i "s/^5$/${NZ}/g" n2to3_cmd.txt;
n2to3 < n2to3_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  genmap ..."
## genmap
echo "${CNAME}_${NZ}z_3D" > genmap_cmd.txt
echo "0.000001" >> genmap_cmd.txt
sed -i "s/2D_//g" genmap_cmd.txt
genmap < genmap_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "done."
# cleanup
rm -f fort.99 gmsh2nek_cmd.txt n2to3_cmd.txt re2torea_cmd.txt genmap_cmd.txt
