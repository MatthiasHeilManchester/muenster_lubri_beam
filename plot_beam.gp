

set xrange [-2:12]
set yrange [-4:10]
set size ratio -1


do for [i=0:1029] { plot "RESLT/beam".i.".dat" u 1:2 w linesp; pause 0.005 }
