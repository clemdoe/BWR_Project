#set nonlogscale y
set title "Flujo de calor entre cada patch"
set ylabel 'Flujo calor [W]'
set xlabel 'Iteration'
#set yrange [-25000:25000]
#set xrange [5000:]
plot "< cat log |grep 'Patch name: solido_to_region0' | cut -d' ' -f6 | tr -d =,= " title 'Solido-DosFluidos' with lines,\
"< cat log |grep 'Patch name: region0_to_solido' | cut -d' ' -f7 | tr -d =,= " title 'DosFluidos-Solido' with lines,\
"< cat log |grep 'Patch name: tubo_to_solido' | cut -d' ' -f7 | tr -d =,= " title 'Primario-Solido' with lines,\
"< cat log |grep 'Patch name: solido_to_tubo' | cut -d' ' -f7 | tr -d =,= " title 'Solido-Primario' with lines

pause 5
reread

