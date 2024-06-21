

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1

if (exist("png")){set terminal png}


if (exist("title")){
  set title title
}

# Validation
do for [i=0:nstep]{
if (exist("png")){set output sprintf('beam%05d.png',i)};
plot "beam".i.".dat" u 1:4 w lines lw 3 lc "black" title "film height step".i."" , "beam".i.".dat" u 1:9 w lines lw 2 lc "green" dt 2 title "exact step".i."";
pause 0.005
}
