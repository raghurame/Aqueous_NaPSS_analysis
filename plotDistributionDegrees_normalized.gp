set term qt font " DejaVu Math TeX Gyre, 14"
set ylabel "Normalized frequency"
set xlabel "{/Symbol q} (degrees)"
set xtics ("0" 0, "30" 10, "60" 20, "90" 30, "120" 40, "150" 50, "180" 60)
plot [:][0:1.4] './degrees.dist.norm' u 1 w l lw 2 title "0x10^{-10} m < r < 4x10^{-10} m", './degrees.dist.norm' u 20 w l lw 2 title "Bulk water"
