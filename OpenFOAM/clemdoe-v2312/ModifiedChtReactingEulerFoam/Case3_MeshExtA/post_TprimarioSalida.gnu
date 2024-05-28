#set nonlogscale y
set title "Temperatura salida fluido monofasico"
set ylabel 'Temperatura'
set xlabel 'Iteration'
set grid y x
#set yrange [0.03125:0.035]
#set xrange [2000:]
plot "< cat log | grep 'weightedAverage(outlet_primario) of T =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida primario' with lines,\
"< cat log | grep 'areaAverage(outlet_secundario) of T.water =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Water' with lines,\
"< cat log | grep 'areaAverage(outlet_secundario) of T.steam =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Steam' with lines

#plot "< cat log | grep 'areaAverage(outlet_secundario) of T.water =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Water' with lines
#plot "< cat log | grep 'areaAverage(outlet_secundario) of T.steam =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Steam' with lines
#plot "< cat log | grep 'areaAverage(outlet_primario) of T =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida primario' with lines,\
#"< cat log | grep 'areaAverage(outlet_secundario) of T.water =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Water' with lines,\
#"< cat log | grep 'areaAverage(outlet_secundario) of T.steam =' | cut -d ' ' -f9 | tr -d =,=" title 'T salida secundario Steam' with lines,\

pause 5
reread


