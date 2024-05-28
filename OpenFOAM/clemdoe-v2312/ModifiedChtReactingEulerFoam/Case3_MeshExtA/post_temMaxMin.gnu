#set nonlogscale y
set title "Temperatura maximas y minimas"
set ylabel 'Temperatura'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log |grep 'max(T.steam)' | cut -d' ' -f7 | tr -d =,=" title 'T.steam.max' with lines,\
"< cat log |grep 'min(T.steam)' | cut -d' ' -f7 | tr -d =,=" title 'T.steam.min' with lines,\
"< cat log |grep 'max(T.water)' | cut -d' ' -f7 | tr -d =,=" title 'T.water.max' with lines,\
"< cat log |grep 'min(T.water)' | cut -d' ' -f7 | tr -d =,=" title 'T.water.min' with lines
pause 5
reread

