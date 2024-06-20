#! /bin/bash


# Movie or on-screen output?
make_movie=1


important_files="muenster_lubri_beam muenster_lubri_beam.cc plot_muenster_lubri_beam.gp validate_plot_muenster.gp plot_muenster.gp"

# Build the thing
make muenster_lubri_beam


# Create run directory
run_dir="Runs"
if [ -e $run_dir ]; then
    echo "Directory "$run_dir" already exists"
    exit 1;
fi

# Move to run dir
mkdir $run_dir
cp $important_files $run_dir
cd $run_dir

full_path_to_run_dir=$PWD


# Do real run (with FSI; try q_fsi = 1.0e-8; 1.0e-7; 5.0e-7)
#===========================================================
echo " " 
echo "----------------------------------------------------------" 
echo " "
reslt_dir=RESLT_FSI
echo "Doing "$reslt_dir
mkdir $reslt_dir
./muenster_lubri_beam --ntstep 500 --t_max 50.0 --q_fsi_target 5.0e-7 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1
if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "

cd $full_path_to_run_dir
# exit 0

# Do real run (zero FSI; pressure kick to check energy conservation; 100 steps up to t=10.0 is OK; 1000 is perfect
#=================================================================================================================
echo " " 
echo "----------------------------------------------------------" 
echo " "
reslt_dir=RESLT_KICK_NO_FSI
echo "Doing "$reslt_dir
mkdir $reslt_dir
./muenster_lubri_beam --ntstep 100 --t_max 10.0 --t_switch_off_kick 2.0 --p_ext_kick 1.0e-5 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1

if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " " 

cd $full_path_to_run_dir
exit 0


# Do real run (zero FSI)
#=======================
reslt_dir=RESLT_NO_FSI
mkdir $reslt_dir
./muenster_lubri_beam --ntstep 500 --t_max 50.0 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1
gnuplot -p -e "nstep=$nstep" ../plot_muenster.gp


cd $full_path_to_run_dir
exit 0


# Do real run (steady only)
#==========================
reslt_dir=RESLT_STEADY
mkdir $reslt_dir
./muenster_lubri_beam --steady_only --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1
gnuplot -p -e "nstep=$nstep" ../plot_muenster.gp

cd $full_path_to_run_dir
exit 0


# Do validation run
#===================
reslt_dir=RESLT_VALIDATION
mkdir $reslt_dir
nstep=100
./muenster_lubri_beam --validate --linearised_validation --ntstep $nstep --t_max 3.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
gnuplot -p -e "nstep=$nstep" ../validate_plot_muenster.gp


cd $full_path_to_run_dir
exit 0






