#! /bin/bash


# Movie or on-screen output?
make_movie=1


important_files="muenster_lubri_beam muenster_lubri_beam.cc plot_muenster_lubri_beam.gp validate_plot_muenster.gp plot_muenster.gp plot_displ_and_energy.gp"

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


######################################################################
#
#  Do stuff
#
######################################################################

# Do real run (zero FSI; pressure kick to check energy conservation; 100 steps up to t=10.0 is OK; 1000 is perfect
#=================================================================================================================
echo " " 
echo "----------------------------------------------------------" 
echo " "
reslt_dir=RESLT_KICK_ONE_WAY_FSI
echo "Doing "$reslt_dir
mkdir $reslt_dir
./muenster_lubri_beam --ntstep 100 --t_max 10.0 --t_switch_off_kick 2.0 --p_ext_kick 1.0e-5 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1

if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; suppress_film=1; title=\"one-way FSI; kick\"" ../plot_muenster.gp
    gnuplot -p -e "nstep=$nstep; png=1; title=\"one-way FSI; kick\"" ../plot_muenster.gp
    #ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    #echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; suppress_film=1; title=\"no FSI; kick\"" ../plot_muenster.gp
fi
gnuplot -p -e "png=1" ../plot_displ_and_energy.gp
echo " " 
echo "----------------------------------------------------------" 
echo " " 

cd $full_path_to_run_dir
exit 0


# Do real run (steady only)
#==========================
reslt_dir=RESLT_STEADY
echo "Doing "$reslt_dir
mkdir $reslt_dir
./muenster_lubri_beam --steady_only --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1

if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; title=\"steady\"" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"steady\"" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "


cd $full_path_to_run_dir
exit 0


# Do real run (zero FSI; pinned film at zero)
#============================================
reslt_dir=RESLT_NO_FSI_LUBRI_PINNED
echo "Doing "$reslt_dir
mkdir $reslt_dir
#./muenster_lubri_beam --ntstep 500 --t_max 50.0 --q_fsi_target 0.0 --nelement 10 --pin_h_lubri_at_zero --dir_name $reslt_dir > $reslt_dir/OUTPUT
./muenster_lubri_beam --ntstep 500 --t_max 1.0 --q_fsi_target 0.0 --linearised_flux --nelement 10 --pin_h_lubri_at_zero --dir_name $reslt_dir > $reslt_dir/OUTPUT

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1

if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; title=\"no FSI; pinned lubri\"" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"no FSI; pinned_lubri\"" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "


cd $full_path_to_run_dir
# exit 0


# Do real run (zero FSI)
#=======================
reslt_dir=RESLT_NO_FSI
echo "Doing "$reslt_dir
mkdir $reslt_dir
./muenster_lubri_beam --ntstep 500 --t_max 50.0 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT
#./muenster_lubri_beam --ntstep 500 --t_max 500.0 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 
#./muenster_lubri_beam --ntstep 5000 --t_max 500.0 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT
#./muenster_lubri_beam --ntstep 500 --t_max 500.0 --nelement 100 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT
#./muenster_lubri_beam --ntstep 500 --t_max 5.0 --q_fsi_target 0.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT
#./muenster_lubri_beam --ntstep 500 --t_max 5.0 --q_fsi_target 0.0 --nelement 51 --dir_name $reslt_dir > $reslt_dir/OUTPUT 




# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1

if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; title=\"no FSI\"" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"no FSI\"" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "


cd $full_path_to_run_dir
# exit 0



# Do validation run with nonlinear flux
#======================================
echo " " 
echo "----------------------------------------------------------" 
echo " "
reslt_dir=RESLT_NONLINEAR_VALIDATION
echo "Doing "$reslt_dir
mkdir $reslt_dir
nstep=100
./muenster_lubri_beam --validate --ntstep $nstep --t_max 3.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1
if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; title=\"nonlinear validation\"" ../validate_plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"nonlinear validation\"" ../validate_plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "

cd $full_path_to_run_dir
# exit 0




# Do validation run with linearised flux
#=======================================
echo " " 
echo "----------------------------------------------------------" 
echo " "
reslt_dir=RESLT_LINEARISED_VALIDATION
echo "Doing "$reslt_dir
mkdir $reslt_dir
nstep=100
./muenster_lubri_beam --validate --linearised_flux --ntstep $nstep --t_max 3.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
nstep=`find . -name 'beam*.dat' | wc -w `
let nstep=$nstep-1
if [ $make_movie -eq 1 ]; then
    gnuplot -p -e "nstep=$nstep; png=1; title=\"linearised validation\"" ../validate_plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"linearised validation\"" ../validate_plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "

cd $full_path_to_run_dir
#exit 0


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
    gnuplot -p -e "nstep=$nstep; png=1; title=\"full FSI\"" ../plot_muenster.gp
    ffmpeg -hide_banner -loglevel error -framerate 5 -pattern_type glob -i 'beam*.png' -c:v ffv1 beam.avi
    echo "Movie in "$PWD"/beam.avi"
else
    gnuplot -p -e "nstep=$nstep; title=\"full FSI\"" ../plot_muenster.gp
fi
echo " " 
echo "----------------------------------------------------------" 
echo " "

cd $full_path_to_run_dir
# exit 0









