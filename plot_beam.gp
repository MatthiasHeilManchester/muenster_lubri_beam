

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1

#set yrange [-0.1:0.1]

#do for [i=0:0] { plot "RESLT/beam".i.".dat" u 1:2 w linesp, "RESLT/beam".i.".dat" u 1:3 w linesp; ; pause 0.005 }
#pause -1



#do for [i=0:3000] { plot "RESLT/beam".i.".dat" u 1:5 w linesp ; pause 0.005 } #, "RESLT/beam".i.".dat" u 1:5 w linesp; ; pause 0.005 }


# FSI
#do for [i=0:7700] { plot "RESLT/beam".i.".dat" u 1:2 w lines lw 2 lc "black" title "beam step".i."" , "RESLT/beam".i.".dat" u 1:3 w lines lw 2 lc "blue" title "film step".i.""; ; pause 0.005 }



#do for [i=31:10029] { plot "RESLT/beam".i.".dat" u 1:4 w linesp, "RESLT/beam".i.".dat" u 1:4 w linesp; ; pause 0.005 }



# Validation
do for [i=0:7700] { plot "RESLT/beam".i.".dat" u 1:4 w lines lw 2 lc "black" title "film height step".i."" , "RESLT/beam".i.".dat" u 1:9 w lines lw 2 lc "blue" title "exact step".i.""; ; pause 0.005 }

# src fct
#do for [i=0:7700] { plot "RESLT/beam".i.".dat" u 1:10 w lines lw 2 lc "black" title "src fct step".i.""; pause 0.005 }
