plot() {
    gnuplot -e "set terminal png; set output 'plot_$1.png'; set xrange [390:800]; plot 'result_$1' using (column(0)*0.125):(column(4)/0.003989422804014327) w l title 'From simulation'"
}

#for i in $(seq 7); do
#    plot $i
#done
plot 1
