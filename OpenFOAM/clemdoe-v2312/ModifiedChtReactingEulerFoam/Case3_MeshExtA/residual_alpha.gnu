#set nonlogscale y
set title "Volumen fraction alpha"
set ylabel 'Residual'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
plot "< cat log |grep 'alpha.steam volume fraction = ' | cut -d' ' -f5 | tr -d =,=" title 'alpha.air volumen fraction' with lines
pause 5
reread

