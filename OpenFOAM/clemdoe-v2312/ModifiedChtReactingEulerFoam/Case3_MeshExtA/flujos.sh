rm -r Flujos
mkdir Flujos
cat log |grep 'Patch name: solido_to_region0' | cut -d' ' -f6 | tr -d =,= > Flujos/SolidoDosFluidos.txt
cat log |grep 'Patch name: region0_to_solido' | cut -d' ' -f7 | tr -d =,= > Flujos/DosFluidosSolido.txt
cat log |grep 'Patch name: anulo_to_solido' | cut -d' ' -f7 | tr -d =,= > Flujos/PrimarioSolido.txt
cat log |grep 'Patch name: solido_to_anulo' | cut -d' ' -f7 | tr -d =,= > Flujos/SolidoPrimario.txt
cat log |grep "^Time =" | cut -d "=" -f2 > Flujos/Time.txt


