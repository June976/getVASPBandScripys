#!/bin/bash

#########################################
### Get Vasp Bands Data, Spin Version ###
#########################################

KPOITS_FILE_NAME=KPOINTS         # KPOINTS file name
OUTCAR_FILE_NAME=OUTCAR		 # OUTCAR file name
FERMI_ENERGY=4.4017              # fermi energy

rec_a_x=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==1{print $4*2.0*3.141592654}')
rec_a_y=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==1{print $5*2.0*3.141592654}')
rec_a_z=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==1{print $6*2.0*3.141592654}')

rec_b_x=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==2{print $4*2.0*3.141592654}')
rec_b_y=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==2{print $5*2.0*3.141592654}')
rec_b_z=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==2{print $6*2.0*3.141592654}')

rec_c_x=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==3{print $4*2.0*3.141592654}')
rec_c_y=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==3{print $5*2.0*3.141592654}')
rec_c_z=$(grep "reciprocal lattice vectors" $OUTCAR_FILE_NAME -A 3|tail -3|awk 'NR==3{print $6*2.0*3.141592654}')

NBANDS=$(awk '/NBANDS/ {print int($NF)}' $OUTCAR_FILE_NAME)
NKPTS=$(awk '/NKPTS/ {print int($4)}' $OUTCAR_FILE_NAME)
NBRCHKPTS=$(cat $KPOITS_FILE_NAME | awk 'NR==2{print int($1)}')
NBRANCH=$(($NKPTS/$NBRCHKPTS))


line_count=$((($NBANDS+3)*$NKPTS))

grep "k-point   1" $OUTCAR_FILE_NAME -A $line_count | head -$line_count > up_bds.txt
grep "k-point   1" $OUTCAR_FILE_NAME -A $(($line_count-1)) | tail -$line_count > dn_bds.txt

function getDistanceOfTwoPoint(){
	### p1_xyz是第一个点坐标，p2_xyz是第二个点的坐标,主要用作获取两个坐标点之间的距离
	dp_x=$(echo $p1_xyz $p2_xyz | awk '{print $4-$1}')
	dp_y=$(echo $p1_xyz $p2_xyz | awk '{print $5-$2}')
	dp_z=$(echo $p1_xyz $p2_xyz | awk '{print $6-$3}')
	dp_xx=$(echo $dp_x $dp_y $dp_z $rec_a_x $rec_b_x $rec_c_x | awk '{print $1*$4+$2*$5+$3*$6}')
	dp_yy=$(echo $dp_x $dp_y $dp_z $rec_a_y $rec_b_y $rec_c_y | awk '{print $1*$4+$2*$5+$3*$6}')
	dp_zz=$(echo $dp_x $dp_y $dp_z $rec_a_z $rec_b_z $rec_c_z | awk '{print $1*$4+$2*$5+$3*$6}')
	distance=$(echo $dp_xx $dp_yy $dp_zz | awk '{print sqrt($1*$1+$2*$2+$3*$3)}')
	echo $distance
}

echo "branch_idx k1_x k1_y k1_z k2_x k2_y k2_z k_label_st k_label_ed pos_x_st pos_x_ed" > klabel.txt
sum_distance=0.0
for i in $(seq 1 $NBRANCH)
do
	nline1=$((i*2-1))
	nline2=$((i*2))
	p1_xyz=$(grep 0\\. $KPOITS_FILE_NAME | awk 'NR==vline{print $1,$2,$3}' vline=$nline1)
	k_label_st=$(grep 0\\. $KPOITS_FILE_NAME | awk 'NR==vline{printf("%s",$NF)}' vline=$nline1)
	p2_xyz=$(grep 0\\. $KPOITS_FILE_NAME | awk 'NR==vline{print $1,$2,$3}' vline=$nline2)
        k_label_ed=$(grep 0\\. $KPOITS_FILE_NAME | awk 'NR==vline{printf("%s",$NF)}' vline=$nline2)
	distance=$(getDistanceOfTwoPoint)
	pos_x_st=$sum_distance
	sum_distance=$(echo $sum_distance $distance | awk '{print $1+$2}')
	echo $i $p1_xyz $p2_xyz $k_label_st $k_label_ed $pos_x_st $sum_distance >> klabel.txt
done

echo "kps up_energy(eV) dn_energy(eV) up_energy_set_zero_fermi dn_energy_set_zero_fermi" > bds_energy.txt
for i in $(seq 1 $NBANDS)
do
	echo "Now deal with Band no $i"
	for j in $(seq 1 $NBRANCH)
	do
		posx_st=$(grep 0 klabel.txt | awk 'NR==vline{print $10}' vline=$j)
		posx_ed=$(grep 0 klabel.txt | awk 'NR==vline{print $11}' vline=$j)
		for k in $(seq 1 $NBRCHKPTS)
		do
			nline=$(($i+2+(($j-1)*$NBRCHKPTS+$k-1)*($NBANDS+3)))
			#echo $nline
			up_energy=$(sed -n "$nline p" up_bds.txt | awk '{print $2}')
			dn_energy=$(sed -n "$nline p" dn_bds.txt | awk '{print $2}')
			posx=$(echo $posx_st $posx_ed $k $NBRCHKPTS | awk '{print $1+($2-$1)*($3-1)/$4}')
			str_one_line=$(echo $posx $up_energy $dn_energy $FERMI_ENERGY | awk '{print $1,$2,$3,$2-$4,$3-$4}')
			echo $str_one_line >> bds_energy.txt
		done	
	done
	echo ' ' >> bds_energy.txt
done

rm up_bds.txt
rm dn_bds.txt
