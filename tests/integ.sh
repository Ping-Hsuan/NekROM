$MOR_DIR/tests/test_template.sh rom_update

ls ../../data/cyl_rect/cyl0.f* > file.list
ls ../../data/cyl_rect/avgcyl0.f* > avg.list
ls bas/bascyl0.f* > bas.list

cp ../../data/cyl_rect/cyl0.f01000 r0.f00001

$SOURCE_ROOT/bin/makenek test
$SOURCE_ROOT/bin/genmap << Z
test
.01
Z

sed -i.bu "s/lb=.*\)/lb=$LB)/g" LMOR
sed -i.bu "s/^.*p177.*\$/$NB p177/g" LMOR

mpiexec -np $NP ./nek5000 | tee test.log.$NP

../$SCR
