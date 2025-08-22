#!/bin/bash


############### INPUT density AND viscosity FILES ###############

#runid=output_test1

#for j in 2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201	202	203	204	205	206	207	208	209	210	211	212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227	228	229	230	231	232	233	234	235	236 ; do

j=147

runid=output_test$j
echo "==================================================================================="
echo "working on RUN ID = $runid..."
echo "==================================================================================="

############### OUTPUT DIRECTORY ###############

outputdir=./MY_OUTPUT
outputdirslash=./MY_OUTPUT/



############### EMPTY DIRECTORIES BEFORE CALCULATION ###############

echo -e "freshening up..."
make realclean
rm -f $outputdir/*
rm -f $outputdir/KERNELS/*
rm -f $runid/*
rm -f $runid/KERNELS/*
rm -f temp.grd
rm -f kernels_grid.xyz
echo -e "done freshening up. \n"



############### POPULATE constants.h file ###############

echo -e "populating constants.h..."

constfile="./const.dat"


echo -e "NOTE, if script freezes here, likely that $constfile does not contain an entry for output test number $j..."

density=`awk -v runid=$runid '{if ($1==runid) print $2}' ${constfile}` # density structure (spherical harmonic, from degree zero up)

visc=`awk -v runid=$runid '{if ($1==runid) print $3}' ${constfile}` # viscosity structure (relative to reference value)

maxdeg=`head -n1 $density | awk '{print $3}' | sed 's/,//g'` # max l for density file (read from file header automatically)

maxdeg_calc=`awk -v runid=$runid '{if ($1==runid) print $4}' ${constfile}` # max spherical harmonic degree FOR DT CALCULATION

mindepth=`awk -v runid=$runid '{if ($1==runid) print $5}' ${constfile}` # min depth in km to use for DT calcs - top_depth in constants.h

kappaindex=`awk -v runid=$runid '{if ($1==runid) print $6}' ${constfile}` # 0 if incompressible (model0, model1), 1 if compressible (sia example)

if [ $kappaindex == "1" ] ; then
	comp="compressible"
else
	comp="incompressible"
fi

rho_w=`awk -v runid=$runid '{if ($1==runid) print $7}' ${constfile}` # air/water loaded

nl=`wc -l $visc | awk '{print $1-1}'` # number of layers, nl+1 = number of radial grid points

visc0=`awk -v runid=$runid '{if ($1==runid) print $8}' ${constfile}`

l_min=`awk -v runid=$runid '{if ($1==runid) print $9}' ${constfile}` # mostly 1, trying 0 to see if it fixes the difference to the TERRA outputs

bgrho=`awk -v runid=$runid '{if ($1==runid) print $10}' ${constfile}` # checking w James Panton which value is appropriate

grav10=`awk -v runid=$runid '{if ($1==runid) print $11}' ${constfile}` # if 1, then g set to =10 everywhere. If not, then g(r) calculated within PM

bdry=`awk -v runid=$runid '{if ($1==runid) print $12}' ${constfile}`

bdry_bot=`awk -v runid=$runid '{if ($1==runid) print $13}' ${constfile}`

echo -e "done populating constants.h. \n"

echo -e "
density = $density
viscosity = $visc
maxdeg (density) = $maxdeg
maxdeg (DT calculation) = $maxdeg_calc
min. depth for DT calc = $mindepth
nl = $nl
kappaindex = $kappaindex, i.e. $comp
rho_w = $rho_w
visc0 = $visc0
l_min = $l_min
bgrho = $bgrho
grav10 = $grav10
bdry = $bdry
bdry_bot = $bdry_bot
"

echo -e "writing constants.h..."
./write_constants_h.sh $maxdeg $nl $maxdeg_calc $mindepth $kappaindex $rho_w $visc0 $l_min $bgrho $grav10 $bdry $bdry_bot
echo -e "done writing constants.h. \n"



######################################

echo -e "running make..."
make 
echo -e "done running make. \n"



######################################

echo -e "running GEOID..."
./GEOID $density $outputdirslash $visc 
echo -e "done running GEOID. \n"



######################################
newoutputdir=./${runid}_LMAX${maxdeg_calc}_RHOW${rho_w}_MIND${mindepth}
newoutputdirnoslash=${runid}_LMAX${maxdeg_calc}_RHOW${rho_w}_MIND${mindepth}

rm -rf $newoutputdir
echo -e "copying output to $newoutputdir..."
mkdir $newoutputdir
cp -r ./MY_OUTPUT/* $newoutputdir
echo -e "done copying. \n"



######################################

# echo -e "plotting..."
# ./plot_quick_dt_output_fixed_ggr.gmt6 $newoutputdirnoslash
# echo -e "done plotting. \n"



############ DONE ###################

done
