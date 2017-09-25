set xtics 90
set xrange [-180:180]
set xlabel "{/Symbol F}"

set ytics 90
set yrange [-180:180]
set ylabel "{/Symbol Y}"

set arrow from -180,0 to 180,0 nohead
set arrow from 0,-180 to 0,180 nohead

unset key
set size square

set title 'Ramachandran Plot'


rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)

plot 'rama.dat' using 2:3:(rgb($4,$5,$6)) linecolor rgb variable pointtype 6 pointsize 1,\
     'rama.dat' using 2:3:(rgb($7,$8,$9)) linecolor rgb variable pointtype 8 pointsize 1,\
     'rama.dat' using 2:3:(rgb($10,$11,$12)) linecolor rgb variable pointtype 4 pointsize 1,\
     'rama.dat' using 2:3:(rgb($13,$14,$15)) linecolor rgb variable pointtype 14 pointsize 1

pause -1
