set term postscript eps enhanced
set output "TotalRuntime.eps"

set yrange[0:2000]
set xlabel "Lattice-Size"
set ylabel "Total Runtime [s]"

plot "./TotalRuntime.txt" u 1:($2/1000000) title "beta = 2","./TotalRuntime.txt" u 1:($3/1000000) title "beta = 8"
