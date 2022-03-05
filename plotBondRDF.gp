set term qt font " DejaVu Math TeX Gyre, 14"
set ylabel "g(r)"
set xlabel "Distance, r (10^{-10} m)"
set xtics nomirror
set ytics nomirror
plot [0:25][:] './bondRDF.output' u 1:2 w lp lw 3 lt 3 title "g_{SO-HO}(r)"
