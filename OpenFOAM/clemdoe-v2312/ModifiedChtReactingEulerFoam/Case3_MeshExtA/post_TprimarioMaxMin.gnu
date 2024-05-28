#set nonlogscale y
set title "Temperatura maximas y minimas - Primario"
set ylabel 'Temperatura'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log |grep 'fieldMinMax minmax-tubo write:'  -A 5| grep 'min(T)' | cut -d' ' -f7 | tr -d =,=" title 'Tmin_primario' with lines,\
"< cat log |grep 'fieldMinMax minmax-tubo write:'  -A 6| grep 'max(T)' | cut -d' ' -f7 | tr -d =,=" title 'Tmax_primario' with lines
pause 5
reread

