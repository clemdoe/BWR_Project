#set nonlogscale y
set title "Fraccion alpha salida"
set ylabel 'alpha'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log | grep 'average(outlet_secundario) of alpha.steam =' | cut -d ' ' -f9 | tr -d =,=" title 'Vapor Salida Secundario' with lines
pause 5
reread
