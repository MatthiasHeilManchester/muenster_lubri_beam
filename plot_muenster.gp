

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1


# beam and free surface
do for [i=0:nstep] { plot "beam".i.".dat" u 1:2 w lines lw 2 lc "black" title "beam step".i."" , "beam".i.".dat" u 1:3 w lines lw 2 lc "blue" title "film step".i.""; ; pause 0.005 }

