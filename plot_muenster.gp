

set xrange [-0.1:1.1]
set yrange [-0.2:0.3]
set size ratio -1



if (exist("title")){
  set title title
}

do for [i=0:nstep]{
 if (exist("png")){if (exist("png")){set terminal png enhanced crop}; set output sprintf('beam%05d.png',i)}

 # Extract time from column 2 of trace file 
 system_command="awk -v step=".i." '\{if ($1==step)\{printf(\"%10.7g\", $2)\}\}' trace_beam.dat" #"echo ".i
 time=system(system_command);

 # Extract pressure from column 3 of trace file 
 system_command="awk -v step=".i." '\{if ($1==step)\{printf(\"%10.7g\", $3)\}\}' trace_beam.dat" 
 pressure=system(system_command);

 # Extract prestress from column 4 of trace file 
 system_command="awk -v step=".i." '\{if ($1==step)\{printf(\"%10.7g\",$4)\}\}' trace_beam.dat" 
 sigma=system(system_command);

 set title "time = ".time."\npre-stress = ".sigma."\n pressure = ".pressure
if ( exist("suppress_film")){
   plot "beam".i.".dat" u 1:2 w lines lw 2 lc "black" title "beam".i.""}
if (!exist("suppress_film")){  
   plot "beam".i.".dat" u 1:2 w lines lw 2 lc "black" title "beam".i."" , "beam".i.".dat" u 1:3 w lines lw 2 lc "blue" title "film".i.""}
 pause 0.005

 if (exist("png")){set terminal pdf enhanced crop; set output sprintf('beam%05d.pdf',i);  replot}
 if (exist("png")){set terminal pdf enhanced crop; set output sprintf('beam_no_padded_number%i.pdf',i);  replot}


}

