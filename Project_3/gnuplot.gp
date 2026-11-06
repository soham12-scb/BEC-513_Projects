nuplot -e "
set terminal pngcairo enhanced font 'Arial,12' size 1000,650;
set output 'vplot_gnuplot.png';
set title 'V-Plot (Gnuplot)';
set xlabel 'Distance from Landmark Center (bp)';
set ylabel 'Fragment Length (bp)';
set xrange [-500:499];
set yrange [50:400];
set view map;
set palette defined (0 'black', 1 'dark-red', 2 'orange', 3 'yellow', 4 'white');
set cblabel 'Tag Density';
plot '< tail -n +2 vplot_matrix.txt | cut -f2-' matrix using (\$1-500):(\$2+50):3 notitle with image;
"
