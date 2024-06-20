

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1

if (exist("png")){set terminal png}

set title GPVAL_PWD
do for [i=0:nstep]{
 if (exist("png")){set output sprintf('beam"%05d".png',i)}
 plot "beam".i.".dat" u 1:2 w lines lw 2 lc "black" title "beam step".i."" , "beam".i.".dat" u 1:3 w lines lw 2 lc "blue" title "film step".i.""
 pause 0.005
}

