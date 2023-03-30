set table "13621060_Hafizh Renanto Akhmad_HW 5.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set xrange[-1:3]; set contour base; set cntrparam levels discrete 0.0; unset surface; set view map; set isosamples 500; splot 2*x*x-2*x*y+2*y*y-1; 
