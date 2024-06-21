

set xrange [-0.1:1.1]
#set yrange [-0.2:0.3]
#set size ratio -1

if (exist("png")){set terminal png}

if (exist("title")){
  set title title
}

do for [i=0:nstep]{
 if (exist("png")){set output sprintf('h_lubri"%05d".png',i)}
 plot "beam".i.".dat" u 1:4 w lines lw 2 lc "blue" title "film step".i.""
 pause -1 # 0.05
}

