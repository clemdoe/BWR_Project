      SUBROUTINE EosLoadIAPWS(fid_lcl)
!
         USE EosIAPWSData
         USE IntrType, ONLY: sdk, sik
         USE Io
!
         IMPLICIT NONE
!
         INTEGER(sik), INTENT(IN) :: fid_lcl        
!
!   This routine was generated by the GenSteamTableRoutines program.
!   This is the driver routine for loading the steam tables from ascii source files.
!   This source file is based on the steam table with the following header:
! tpfh2ogibb version 1.0, tables of thermodynamic properties of light wate
! generated on 28-Nov-07 at 11:17:19 by stgh2onew  1.0 (1/Mar/01)
!
         INTEGER(sik) :: ntot, nuse, i
!
         CALL EosLoadIAPWSScalars(ttrip, ptrip, vtrip, &
              tcrit, pcrit, vcrit, tmin, pmin, tmax, pmax, &
              ptable1, ltable1, ntable1, stable1, ptable2, ltable2, ntable2, stable2, &
              ptable3, ltable3, ntable3, stable3, ptable4, ltable4, ntable4, stable4, &
              ptable5, ltable5, ntable5, stable5, ptable6, ltable6, ntable6, stable6, &
              ptable7, ltable7, ntable7, stable7, ptable8, ltable8, ntable8, stable8, &
              ptable9, ltable9, ntable9, stable9, ptable10, ltable10, ntable10, stable10, &
              ntemp, npres, nsubcrittemp, nsupcrittemp, nsattemp, ofirstsattemp, &
              olastsattemp, otemptrip1, otempcrit1, nsubcritpres, nsupcritpres, &
              nsatpres, ofirstsatpres, olastsatpres, oprestrip2, oprescrit2, &
              ofirstsupcritpres, olastsupcritpres, ntot, newrecord)
!
         DO i = 1, SIZE(outunt)
!            WRITE (outunt(i), '(a)') ' Default steam tables with the following header will be used:'
!            WRITE (outunt(i), '(1x, a/1x, a)') newrecord(1:2)
            WRITE (fid_lcl, '(a)') ' Default steam tables with the following header will be used:'
            WRITE (fid_lcl, '(1x, a/1x, a)') newrecord(1:2)
         END DO
!
         CALL AllocIAPWSTables(ntot, nuse)
         IF (nuse == -1) THEN
             WRITE(*,*) '*EosLoadIAPWS* error allocating tables'
            !CALL error(1, ' *EosLoadIAPWS* Error allocating tables.')
         END IF
!
         CALL EosLoadIAPWSTable1
         CALL EosLoadIAPWSTable2
         CALL EosLoadIAPWSTable3
         CALL EosLoadIAPWSTable4
         CALL EosLoadIAPWSTable5
         CALL EosLoadIAPWSTable6
         CALL EosLoadIAPWSTable7
         CALL EosLoadIAPWSTable8
         CALL EosLoadIAPWSTable9
         CALL EosLoadIAPWSTable10
         CALL EosLoadIAPWSTable11
         CALL EosLoadIAPWSTable12
         CALL EosLoadIAPWSTable13
         CALL EosLoadIAPWSTable14
         CALL EosLoadIAPWSTable15
         CALL EosLoadIAPWSTable16
         CALL EosLoadIAPWSTable17
         CALL EosLoadIAPWSTable18
         CALL EosLoadIAPWSTable19
         CALL EosLoadIAPWSTable20
         CALL EosLoadIAPWSTable21
         CALL EosLoadIAPWSTable22
         CALL EosLoadIAPWSTable23
         CALL EosLoadIAPWSTable24
         CALL EosLoadIAPWSTable25
         CALL EosLoadIAPWSTable26
         CALL EosLoadIAPWSTable27
         CALL EosLoadIAPWSTable28
         CALL EosLoadIAPWSTable29
         CALL EosLoadIAPWSTable30
         CALL EosLoadIAPWSTable31
         CALL EosLoadIAPWSTable32
         CALL EosLoadIAPWSTable33
         CALL EosLoadIAPWSTable34
         CALL EosLoadIAPWSTable35
         CALL EosLoadIAPWSTable36
         CALL EosLoadIAPWSTable37
         CALL EosLoadIAPWSTable38
         CALL EosLoadIAPWSTable39
         CALL EosLoadIAPWSTable40
         CALL EosLoadIAPWSTable41
         CALL EosLoadIAPWSTable42
         CALL EosLoadIAPWSTable43
         CALL EosLoadIAPWSTable44
         CALL EosLoadIAPWSTable45
         CALL EosLoadIAPWSTable46
         CALL EosLoadIAPWSTable47
         CALL EosLoadIAPWSTable48
         CALL EosLoadIAPWSTable49
         CALL EosLoadIAPWSTable50
!
      END SUBROUTINE EosLoadIAPWS
