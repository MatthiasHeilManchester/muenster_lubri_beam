

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1

#set yrange [-0.1:0.1]

#do for [i=0:0] { plot "RESLT/beam".i.".dat" u 1:2 w linesp, "RESLT/beam".i.".dat" u 1:3 w linesp; ; pause 0.005 }
#pause -1



#do for [i=0:3000] { plot "RESLT/beam".i.".dat" u 1:5 w linesp ; pause 0.005 } #, "RESLT/beam".i.".dat" u 1:5 w linesp; ; pause 0.005 }


do for [i=0:3000] { plot "RESLT/beam".i.".dat" u 1:2 w linesp, "RESLT/beam".i.".dat" u 1:3 w linesp; ; pause 0.005 }

#do for [i=31:10029] { plot "RESLT/beam".i.".dat" u 1:4 w linesp, "RESLT/beam".i.".dat" u 1:4 w linesp; ; pause 0.005 }

