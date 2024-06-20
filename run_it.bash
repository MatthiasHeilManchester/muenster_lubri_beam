#! /bin/bash


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


exit 0;


# Do validation run
#===================
reslt_dir=RESLT_VALIDATION
mkdir $reslt_dir
nstep=100
./muenster_lubri_beam --validate --linearised_validation --ntstep $nstep --t_max 3.0 --dir_name $reslt_dir > $reslt_dir/OUTPUT 

# Gnuplot the the thing
cd $reslt_dir
gnuplot -p -e "nstep=$nstep" ../validate_plot_muenster.gp









