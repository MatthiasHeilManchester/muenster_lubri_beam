

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1

# Validation
do for [i=0:nstep] { plot "beam".i.".dat" u 1:4 w lines lw 2 lc "black" title "film height step".i."" , "beam".i.".dat" u 1:9 w lines lw 2 lc "blue" title "exact step".i.""; ; pause 0.005 }
