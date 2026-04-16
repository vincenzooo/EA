plot 'Shell1_imd.txt'  using 1:($2**2*38.5817083793337) w l,\
'Effective_area_total.dat' using 1:2 w p

set grid

set terminal png
set output 'Shell1_verifica.png'
replot
set terminal window
set output
