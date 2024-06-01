

set xrange [-2:12]
set yrange [-4:10]
set size ratio -1


#do for [i=0:0] { plot "RESLT/beam".i.".dat" u 1:2 w linesp, "RESLT/beam".i.".dat" u 1:3 w linesp; ; pause 0.005 }
#pause -1


#do for [i=0:1029] { plot "RESLT/beam".i.".dat" u 1:2 w linesp, "RESLT/beam".i.".dat" u 1:3 w linesp; ; pause 0.005 }

do for [i=31:10029] { plot "RESLT/beam".i.".dat" u 1:4 w linesp, "RESLT/beam".i.".dat" u 1:4 w linesp; ; pause 0.005 }

