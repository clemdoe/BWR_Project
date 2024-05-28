#set nonlogscale y
set title "Volumen fraction alpha"
set ylabel 'Residual'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log |grep 'Tf.steamAndWater' | cut -d' ' -f4 | tr -d =,=" title 'T.SteamAndWater_min' with lines,\
"< cat log |grep 'Tf.steamAndWater' | cut -d' ' -f7 | tr -d =,=" title 'T.SteamAndWater_averg' with lines,\
"< cat log |grep 'Tf.steamAndWater' | cut -d' ' -f10 | tr -d =,=" title 'T.SteamAndWater_max' with lines
pause 5
reread

