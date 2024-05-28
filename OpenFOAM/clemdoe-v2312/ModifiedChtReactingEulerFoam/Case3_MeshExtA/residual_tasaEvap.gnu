#set nonlogscale y
set title "Tasas y flujos de calor"
set ylabel 'Flujo calor'
set xlabel 'Iteration'
#set yrange [0.03125:0.035]
set grid y x
plot "< cat log |grep 'iDmdt.steamAndWater:' | cut -d' ' -f4 | tr -d =,=" title 'iDmdt.min' with lines,\
"< cat log |grep 'iDmdt.steamAndWater:' | cut -d' ' -f7 | tr -d =,=" title 'iDmdt.medio' with lines,\
"< cat log |grep 'iDmdt.steamAndWater:' | cut -d' ' -f10 | tr -d =,=" title 'iDmdt.max' with lines,\
"< cat log |grep 'iDmdt.steamAndWater:' | cut -d' ' -f13 | tr -d =,=" title 'iDmdt.integral' with lines,\
 "< cat log |grep 'wDmdt.steamAndWater:' | cut -d' ' -f4 | tr -d =,=" title 'wDmdt.min' with lines,\
"< cat log |grep 'wDmdt.steamAndWater:' | cut -d' ' -f7 | tr -d =,=" title 'wDmdt.medio' with lines,\
"< cat log |grep 'wDmdt.steamAndWater:' | cut -d' ' -f10 | tr -d =,=" title 'wDmdt.max' with lines,\
"< cat log |grep 'wDmdt.steamAndWater:' | cut -d' ' -f13 | tr -d =,=" title 'wDmdt.integral' with lines

pause 5
reread

