#set nonlogscale y
set title "Temperatura maximas y minimas"
set ylabel 'Temperatura'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log |grep 'max(mag(U.steam))' | cut -d' ' -f7 | tr -d =,=" title 'U.steam.max' with lines,\
"< cat log |grep 'min(mag(U.steam))' | cut -d' ' -f7 | tr -d =,=" title 'U.steam.min' with lines,\
"< cat log |grep 'min(mag(U.water))' | cut -d' ' -f7 | tr -d =,=" title 'U.water.max' with lines,\
"< cat log |grep 'min(mag(U.water))' | cut -d' ' -f7 | tr -d =,=" title 'U.water.min' with lines
pause 5
reread
