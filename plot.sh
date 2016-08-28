plot() {
    gnuplot -e "set terminal png; set output 'plot_$i.png'; plot 'result_$i' using 4 w l;"
}

for i in $(seq 7); do
    plot $i
done
