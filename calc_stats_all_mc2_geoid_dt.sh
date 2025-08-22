#!/bin/bash

#echo "MODEL GEOID_SCORE RL RL10 L2_NORM L2_NORM_LATFAC CMEAN CSTD RMSP RMSPLOG RMSP10 RMSPLOG10 MC2NAME" > geoid_scores_all_mc2_new.dat

echo "MODEL SIA_GEOID_SCORE CHI_P CHI OVERLINE_RL R M2NAME" > geoid_scores_mc2_4dp.dat
echo "MODEL CHI_P CHI OVERLINE_RL R MC2NAME"  > dt_scores_mc2_4dp.dat


#for runcode in 1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241 ; do

#for runcode in 239 240	241 242 243 244; do

for runcode in 147 ; do

mc2name=`awk -v var=$runcode '{if ($1=="output_test"var) print $2}' < inputs_ggr_all_mc2_runs.dat | sed 's^./INPUT2/DENSITY_INPUTS/^^g' | sed -n  's/_density.*//p'`

echo ""
echo "======== working on" $runcode ":" $mc2name "========="


outdirfull=./output_test${runcode}_LMAX50_RHOW1.03E+03_MIND0   #include all depths

# reference grids and spherical harmonic files
refdtgrd=$outdirfull/GRD_REF_dyntopography
refdtsph=$outdirfull/SPH_REF_dyntopography

refgeogrd=$outdirfull/GRD_REF_geoid
refgeosph=$outdirfull/SPH_REF_geoid

# calculated grids and spherical harmonic files
surftopo=$outdirfull/GRD_surftopo
surftoposph=$outdirfull/SPH_surftopo

geoidgrd=$outdirfull/GRD_geoid
geoidsph=$outdirfull/SPH_geoid


kno=`ls $outdirfull/KERNELS/surftopo* | wc -l | awk '{$1=$1};1'`
echo "LMAX = " $kno


#############################################################################################
# COMPARE DYNAMIC TOPOGRAPHIES
#############################################################################################
echo "--------------- DYNAMIC TOPO STATS ---------------"

# calculate chi (rms) misfit between reference and calculated grids
	f1=$surftopo
	f2=$refdtgrd
	awk '{print cos($2*2*3.14159/360)*(0.5*3.14159)}' $f1 > factors.dat
	awk '{print $3}' $f1 > temp1.dat
	awk '{print $3}' $f2 > temp2.dat
	paste temp1.dat temp2.dat | awk '{print (($1-$2)**2)**0.5}' > temp3.dat
	paste temp3.dat factors.dat | awk '{print $1*$2}' > temp4.dat
	len=`wc -l temp3.dat | awk '{print $1}'`
	l2norm=`awk '{sum+=$1;} END{print sum;}' temp3.dat | awk -v len=$len '{print $1/len}'`
	l2norm_factored=`awk '{sum+=$1;} END{print sum;}' temp4.dat | awk -v len=$len '{print $1/len}'`
	
	l2norm_factored=`echo $l2norm_factored | awk '{printf "%.0f\n", $1}'`


# COMPARE power spectra of calculated surface topography, reference DT and Kaula's rule

# reference power spectrum...
ghead -n-1 $refdtsph | awk '{if (NR>4) print $2, $3, $4, $5}'  > l_m_c1_c2.temp  #GGR adjustment for mac (delete final row command revised)
lmax=`gmt gmtinfo -C l_m_c1_c2.temp | awk '{print $2}'`

	awk '{print $1, $2, $3}' l_m_c1_c2.temp > junk1
	awk '{if ($2>0) print $1, -$2, $4}' l_m_c1_c2.temp > junk2
	cat junk2 >> junk1
	sort -n -k1 -k2 junk1 > junk3
	echo "0 0 0" > junk4 ### ONLY FOR PROPMAT: NO DEGREE 0 STRUCTURE POSSIBLE
	cat junk3 >> junk4
	cp junk4 l_m_coeffs.temp
	rm junk*
	echo -n "" > l_power_ref_dt.temp
	i=0
	while [ $i -le $lmax ] ; do
		awk -v i=$i '{if ($1==i) print $0}' l_m_coeffs.temp > lmc.temp
		power=`awk '{print $3**2.}' lmc.temp | awk '{sum+=$1;} END{printf "%.8f\n", sum;}'`
		echo $i $power >> l_power_ref_dt.temp
		i=$[$i+1]
	done


# calculated power spectrum...
ghead -n-1 $surftoposph | awk '{if (NR>4) print $2, $3, $4, $5}'  > l_m_c1_c2.temp  #GGR adjustment for mac (delete final row)
lmax=`gmt gmtinfo -C l_m_c1_c2.temp | awk '{print $2}'`

	awk '{print $1, $2, $3}' l_m_c1_c2.temp > junk1
	awk '{if ($2>0) print $1, -$2, $4}' l_m_c1_c2.temp > junk2
	cat junk2 >> junk1
	sort -n -k1 -k2 junk1 > junk3
	echo "0 0 0" > junk4 ### ONLY FOR PROPMAT: NO DEGREE 0 STRUCTURE POSSIBLE
	cat junk3 >> junk4
	cp junk4 l_m_coeffs.temp
	rm junk*
	echo -n "" > l_power.temp
	i=0
	while [ $i -le $lmax ] ; do
		awk -v i=$i '{if ($1==i) print $0}' l_m_coeffs.temp > lmc.temp
		power=`awk '{print $3**2.}' lmc.temp | awk '{sum+=$1;} END{printf "%.8f\n", sum;}'`
		echo $i $power >> l_power.temp
		i=$[$i+1]
	done

	
# calculate correlation coefficients, r_l, between reference dynamic topo and calculated surface topo
	echo "0 0 0 0" > t1.temp
	echo "0 0 0 0" > t2.temp
		
	ghead -n-1 $refdtsph | awk '{if (NR>4) print $2, $3, $4, $5}' >> t1.temp # ref
	ghead -n-1 $surftoposph | awk '{if (NR>4) print $2, $3, $4, $5}' >> t2.temp # pred

	python3.11 forte_calc_corr_ggr.py

	gmt math corr1.temp MEAN = m.temp
	gmt math corr1.temp STD = s.temp
	cmean=`gmt gmtinfo -C m.temp | awk '{print $2}'` # conor mean r_l
	cstd=`gmt gmtinfo -C s.temp | awk '{print $2}'` # conor std r_l
	rltot=`awk '{print $0}' < rltot.temp`
	

# kaula's rule
	G="6.67e-11" # m3 kg-1 s-2
	M="5.97e24" # kg
	R="6371000" # m

	for Z in 9 12 15 20 30 40 50 ; do
		Z_adj=`echo $Z | awk '{print $1*1e-8}'` # SI
#		Z_adj=$Z
		k1=`echo $G $M | awk '{print ($1*$2)}'`
		k2=`echo $Z_adj $R | awk '{print ($1*($2**2.))}'`
		k3=`echo $k1 $k2 | awk '{print ($1/$2)**2.}'`
		echo -n "" > kl_Z${Z}.temp
		for i in $(seq 2 1 50) ; do
			kl1=`echo $i | awk '{print 2/$1}'`
			kl2=`echo $i | awk '{print 3/($1**2.)}'`
			kl3=`echo $i | awk '{print 1/($1**4.)}'`
			klfull=`echo $kl1 $kl2 $kl3 | awk '{print $1 - $2 + $3}'`
			echo $i $k3 $klfull | awk '{print $1, ($2*$3)/1e10}' >> kl_Z${Z}.temp # m2
		done
	done
	cat kl_Z9.temp | awk '{print $1, $2/1e6}' > t1.temp
	tac kl_Z15.temp | awk '{print $1, $2/1e6}' >> t1.temp


# 	Holdt 2022 power spectrum
loaded=water

# calculate total pspec misfit (chi_p)
	cp l_power.temp f1
	cat kl_Z12.temp | awk -v lmax=$kno '{if ($1<=lmax) print $1, $2/1e6}' > f2

	if [ $loaded == "water" ] ; then
		#### SEE ipython notebook in from_megan ####
		infiletemp="./from_megan/outputfolder/DT_b_1.00.txt"
	elif [ $loaded == "air" ] ; then
		infiletemp="./from_megan/outputfolder_air/DT_b_1.00.txt"
	fi
	awk -v lmax=$kno '{if ($1<=lmax) print $1, $3}' $infiletemp > f3

	## Kaula misfit
	awk '{if ($1>1) print $01, $2/1e6}' f1 > f1a
	n=`wc -l f1a | awk '{print $1}'`
	rmsk=`paste f1a f2 | awk '{print $2, $4}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`
	rmsklog=`paste f1a f2 | awk '{print log($2)/log(10), log($4)/log(10)}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`

	## Megan misfit 
	awk '{if ($1>0) print $0}' f1 > f1b
	n=`wc -l f1b | awk '{print $1}'`
	rmsm=`paste f1b f3 | awk '{print $2, $4}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`
	rmsmlog=`paste f1b f3 | awk '{print log($2)/log(10), log($4)/log(10)}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`
	
	rmstot=`echo $rmsk $rmsm | awk '{printf "%.1f\n", ($1+$2)}'`
	rmstotlog=`echo $rmsklog $rmsmlog | awk '{printf "%.4f\n", ($1+$2)}'`
	
# trimming to 2 d.p. 	
	cmean=`echo $cmean | awk '{printf "%.4f\n", $1}'`
	rltot=`echo $rltot | awk '{printf "%.4f\n", $1}'`

echo "MODEL CHI_P CHI OVERLINE_RL R MC2NAME" 
echo $runcode $rmstotlog $l2norm_factored $cmean $rltot $mc2name 
echo $runcode $rmstotlog $l2norm_factored $cmean $rltot $mc2name >> dt_scores_mc2_4dp.dat 

#######################################################################################################################################
# COMPARE GEOIDS
#######################################################################################################################################

echo "-------------------- GEOID STATS --------------------"


# SIA'S GEOID SCORE
	infile="$outdirfull/Summary"
	sia_geoid_score=`cat $infile | grep "GeoC" | awk '{print $2}' | sed 's/|//g'`

# calculate chi (rms) misfit between reference and calculated grids
	f1=$geoidgrd
	f2=$refgeogrd
	awk '{print cos($2*2*3.14159/360)*(0.5*3.14159)}' $f1 > factors.dat
	awk '{print $3}' $f1 > temp1.dat
	awk '{print $3}' $f2 > temp2.dat
	paste temp1.dat temp2.dat | awk '{print (($1-$2)**2)**0.5}' > temp3.dat
	paste temp3.dat factors.dat | awk '{print $1*$2}' > temp4.dat
	len=`wc -l temp3.dat | awk '{print $1}'`
	l2norm=`awk '{sum+=$1;} END{print sum;}' temp3.dat | awk -v len=$len '{print $1/len}'`
	l2norm_factored=`awk '{sum+=$1;} END{print sum;}' temp4.dat | awk -v len=$len '{print $1/len}'`


# COMPARE CALCULATED POWER SPECTRUM 
# reference power spectrum...	
ghead -n-1 $refgeosph | awk '{if (NR>4) print $2, $3, $4, $5}'  > l_m_c1_c2.temp  #GGR adjustment for mac (delete final row)
lmax=`gmt gmtinfo -C l_m_c1_c2.temp | awk '{print $2}'`
	awk '{print $1, $2, $3}' l_m_c1_c2.temp > junk1
	awk '{if ($2>0) print $1, -$2, $4}' l_m_c1_c2.temp > junk2
	cat junk2 >> junk1
	sort -n -k1 -k2 junk1 > junk3
	echo "0 0 0" > junk4 ### ONLY FOR PROPMAT: NO DEGREE 0 STRUCTURE POSSIBLE
	cat junk3 >> junk4
	cp junk4 l_m_coeffs.temp
	rm junk*
	
	echo -n "" > l_power_ref_geo.temp
	i=0
	while [ $i -le $lmax ] ; do
		awk -v i=$i '{if ($1==i) print $0}' l_m_coeffs.temp > lmc.temp
		power=`awk '{print $3**2.}' lmc.temp | awk '{sum+=$1;} END{printf "%.8f\n", sum;}'`
		echo $i $power >> l_power_ref_geo.temp
		i=$[$i+1]
	done


# calculated power spectrum...
ghead -n-1 $geoidsph | awk '{if (NR>4) print $2, $3, $4, $5}'  > l_m_c1_c2.temp  #GGR adjustment for mac (delete final row)
lmax=`gmt gmtinfo -C l_m_c1_c2.temp | awk '{print $2}'`
	awk '{print $1, $2, $3}' l_m_c1_c2.temp > junk1
	awk '{if ($2>0) print $1, -$2, $4}' l_m_c1_c2.temp > junk2
	cat junk2 >> junk1
	sort -n -k1 -k2 junk1 > junk3
	echo "0 0 0" > junk4 ### ONLY FOR PROPMAT: NO DEGREE 0 STRUCTURE POSSIBLE
	cat junk3 >> junk4
	cp junk4 l_m_coeffs.temp
	rm junk*

	echo -n "" > l_power.temp
	i=0
	while [ $i -le $lmax ] ; do
		awk -v i=$i '{if ($1==i) print $0}' l_m_coeffs.temp > lmc.temp
		power=`awk '{print $3**2.}' lmc.temp | awk '{sum+=$1;} END{printf "%.8f\n", sum;}'`
		echo $i $power >> l_power.temp
		i=$[$i+1]
	done
	
# calculate correlation coefficients, r_l, between reference and calculated geoid
	echo "0 0 0 0" > t1.temp
	echo "0 0 0 0" > t2.temp
		
	ghead -n-1 $refgeosph | awk '{if (NR>4) print $2, $3, $4, $5}' >> t1.temp # ref
	ghead -n-1 $geoidsph | awk '{if (NR>4) print $2, $3, $4, $5}' >> t2.temp # pred

	python3.11 forte_calc_corr_ggr.py

	gmt math corr1.temp MEAN = m.temp
	gmt math corr1.temp STD = s.temp
	cmean=`gmt gmtinfo -C m.temp | awk '{print $2}'` # conor mean r_l
	cstd=`gmt gmtinfo -C s.temp | awk '{print $2}'` # conor std r_l
	rltot=`awk '{print $0}' < rltot.temp`

	# calculate total pspec misfit, chi_p
	cat l_power.temp | awk '{print $1, $2/1e6}' >  f1
	
	# calc_pow_rms
	cat l_power_ref_geo.temp | awk '{print $1, $2/1e6}' > f2

	## geoid misfit
	awk '{if ($1>1) print $0}' f1 > f1a
	awk '{if ($1>1) print $0}' f2 > f2a
	n=`wc -l f1a | awk '{print $1}'`
	rmsg=`paste f1a f2a | awk '{print $2, $4}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`
	rmsglog=`paste f1a f2a | awk '{print log($2)/log(10), log($4)/log(10)}' | awk '{print ($1-$2)**2.}' | awk '{sum+=$1;}END{print sum;}' | awk -v n=$n '{print ($1/n)**0.5}'`
	
	rmsg=`echo $rmsg | awk '{printf "%.4f\n", $1}'`
	rmsglog=`echo $rmsglog | awk '{printf "%.4f\n", $1}'`

# trimming to 2 d.p. 
	l2norm=`echo $l2norm | awk '{printf "%.4f\n", $1}'`
	l2norm_factored=`echo $l2norm_factored | awk '{printf "%.4f\n", $1}'`
	cmean=`echo $cmean | awk '{printf "%.4f\n", $1}'`
	cstd=`echo $cstd | awk '{printf "%.4f\n", $1}'`
	rltot=`echo $rltot | awk '{printf "%.4f\n", $1}'`
	rl10=`awk '{printf "%.4f", $1}' < rl10.temp`  #rl for all degrees up to l = 10.

# output
# echo "MODEL GEOID_SCORE RL RL10 L2_NORM L2_NORM_LATFAC CMEAN CSTD RMSG RMSGLOG MC2NAME" 
# echo $runcode $sia_geoid_score $rltot $rl10 $l2norm $l2norm_factored $cmean $cstd $rmsg $rmsglog $mc2name # NB stats is 3 cols I think
# echo $runcode $sia_geoid_score $rltot $rl10 $l2norm $l2norm_factored $cmean $cstd $rmsg $rmsglog $mc2name >> geoid_scores_mc2_4dp.dat 

echo "MODEL SIA_GEOID_SCORE CHI_P CHI OVERLINE_RL R M2NAME" 
echo $runcode $sia_geoid_score $rmsglog $l2norm_factored $cmean $rltot $mc2name
echo $runcode $sia_geoid_score $rmsglog $l2norm_factored $cmean $rltot $mc2name >> geoid_scores_mc2_4dp.dat 


done

