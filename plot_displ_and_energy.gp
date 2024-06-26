# Offset in the order left, right, top, bottom
set offset graph 0.0, graph 0.0, graph 0.4, graph 0.05

if (exist("png")){set terminal pdf enhanced crop; set output "displ_and_energy.pdf"}


set xlabel "time"
set ylabel "pressure"
set ytics autofreq
set y2label "centre displacement"
set y2tics autofreq

set multiplot layout 2,1 columns


plot "trace_beam.dat" u 2:3 w lines lw 2 t "pressure", "trace_beam.dat" u 2:5 w lines lw 2 t "centre displacement" axis x1y2

set ylabel "energy"
set ytics autofreq
unset y2label
unset y2tics

plot "trace_beam.dat" u 2:7 w lines lw 2 t "kinetic energy", "trace_beam.dat" u 2:8 w lines lw 2 t "strain energy" , "trace_beam.dat" u 2:($7+$8) w lines lw 2 t "total energy"





unset multiplot