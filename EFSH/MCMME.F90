! THE MODULE OF CALCULATION OF MOLECULAR MATRIX ELEMENTS VER 1.0 11.2005
! бяе опюбю гюыхыемш. йнлоюмхъ LabComputerSistemTecnologyKasprzhitsky
! VER 1.0 NEW  11,2005 цнд

! лндскэ пюяверю лнкейскъпмшу люрпхвмшу щкелемрнб бепяхъ 1.0    
	

module mcmme
implicit none
	
      

contains


      
	
! ондопнцпюллю пюяверю люрпхвмнцн щкелемрю ндмнщкейрпнммни щмепцхх 
! нохяюмхе оюпюлерпнб ондопнцпюллш
! L1-нпахрюкэмши лнлемр оепбнцн щкейрпнмю
! M1-опнейжхъ нпахрюкэмнцн лнлемрю оепбнцн щкейрпнмю
! L2-нпахрюкэмши лнлемр брнпнцн щкейрпнмю
! M2-опнейжхъ нпахрюкэмнцн лнлемрю брнпнцн щкейрпнмю 
! N-вхякн ъдеп
! Z0-гюпъд ъдпю б мювюке йннпдхмюр
! Z(N)-люяяхб гюпъднб ъдеп
! COORRR(N,2)-люяяхб пюдхюкэмшу йнндхмюр 
!  COORRR(N,1)-пюдхюкэмюъ йннпдхмюрю ъдпю (R-пюдхся)
!  COORRR(N,2)-пюдхюкэмюъ йннпдхмюрю ъдпю б мнбни оепелеммни (RO-пюдхся)
! COOAAA(N,2)-люяяхб йннпдхмюр ъдеп
!  COOAAA(N,1)-брнпюъ йннпдхмюрю -Reta-сцнк (0=<Reta<=PI)
!  COOAAA(N,2)-рперэъ йннпдхмюрю -RFu-сцнк (0=<RFu<=2*PI)
! RO0-мювюкэмне гмювемхе хмрепбюкю он RO
! H-ьюц
! Npoint-вхякн рнвей
! R(Npoint)-люяяхб гмювемхи юпцслемрю
! Ffun1(Npoint)-оепбюъ бнкмнбюъ тсмйжхъ F-рхою
! Ffun2(Npoint)-брнпюъ бнкмнбюъ тсмйжхъ F-рхою
! RO1(Npoint)-люяяхб оепеяверю тсмйжхи (ОЕПБЮЪ ОПНХГБНДМНЪ МНБНИ ОЕПЕЛЕММНИ)
! RO2(Npoint)-люяяхб оепеяверю тсмйжхи (БРНПЮЪ ОПНХГБНДМНЪ МНБНИ ОЕПЕЛЕММНИ)
! RO3(Npoint)-люяяхб оепеяверю тсмйжхи (РПЕРЭЪ ОПНХГБНДМНЪ МНБНИ ОЕПЕЛЕММНИ)
! EnergyKIN-йхмерхвеяйюъ вюярэ щмепцхх
! EnergyPOT-онремжхюкэмюъ вюярэ щмепцхх
real(8) function CMME_ONE_ELECTRONIC_INTEGRAL(L1,M1,L2,M2,N,Z0,Z,COORRR,COOAAA,RO0,H,Npoint,R,Ffun1,Ffun2,RO1,RO2,RO3,EnergyKIN,EnergyPOT)
 use mcmri,only:CMRI_NABLA,CMRI_KEN,CMRI_CENT,CMRI_CENT_R
 use mcmav,only:CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD
 implicit none
 integer::L1,M1,L2,M2,N,Npoint
 real(8)::Z0,RO0,H,EnergyKIN,EnergyPOT
 real(8),dimension(:)::Z,R,Ffun1,Ffun2,RO1,RO2,RO3
 real(8),dimension(:,:)::COORRR,COOAAA
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::ierr,IKNumbre,ISCMME,IKK,ISA
 real(8)::ENERGYDIAGONAL,ENABLA,EKIN,ECENT
 real(8)::ENERGYCENTZA,ESUM,EZSUM
 real(8),allocatable,dimension(:,:)::RintK
 real(8),allocatable,dimension(:,:,:)::RcoffACF




 ! щрюо 1. нясыеярбкъел пюявер дхюцнмюкэмни вюярх люрпхвмнцн щкелемрю
 ENERGYDIAGONAL=0.D0
 ENABLA=0.D0
 EKIN=0.D0
 ECENT=0.D0
 IF(L1.EQ.L2.AND.M1.EQ.M2) THEN
   ! люрпхвмши щкелемр брнпни опнхгбндмни 
   ENABLA=-0.5D0*CMRI_NABLA(Npoint,RO0,H,R,Ffun1,Ffun2,RO1,RO2,RO3)
   ! люрпхвмши щкелемр жемрпнаефмнцн вкемю
   EKIN=FLOAT(L1*(L1+1))*0.5D0*CMRI_KEN(RO0,H,Npoint,R,Ffun1,Ffun2,RO1)
   ! люрпхвмши щкелемр ноепюрнпю бгюхлндеиярбхъ щкейрпнмю я ъдпнл мюундъыхляъ б мювюке йннпдхмю
   IF(Z0.NE.0.D0) THEN
      ECENT=-Z0*CMRI_CENT(RO0,H,Npoint,R,Ffun1,Ffun2,RO1) 
     ELSE
	  ECENT=0.D0 
   ENDIF  
   ! пегскэрюр
   ENERGYDIAGONAL=ENABLA+EKIN+ECENT
 ENDIF
 
 ! щрюо 2. нясыеярбкъел пюявер щмепцхх бгюхлндеиярбхъ щкейрпнмю я ъдпюлх ме мюундъыхлхяъ б мювюке йннпдхмюр
 ENERGYCENTZA=0.D0
 IF(N.NE.0) THEN
    ! СЯРЮМЮБКХБЮЕЛ ПЮГЛЕП ЛЮЯЯХБЮ
    IKNumbre=0
    DO ISCMME=IABS(L1-L2),L1+L2,2
       IKNumbre=IKNumbre+1
    ENDDO
    !бшдекъел оюлърэ дкъ люяяхбнб
    allocate(RintK(N,IKNumbre),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CMME_ONE_ELECTRONIC_INTEGRAL'
      write(*,*) 'MEMORY ON THE FILE "RintK" IS NOT SELECTED'
      stop 
    endif 
    allocate(RcoffACF(N,IKNumbre,2),stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CMME_ONE_ELECTRONIC_INTEGRAL'
      write(*,*) 'MEMORY ON THE FILE "RcoffACF" IS NOT SELECTED'
      stop 
    endif 
    ! гюмскъел оепед хяонкэгнбюмхел
    RintK=0.D0
    RcoffACF=0.D0
    ! пюявхршбюел йнщттхжхемрш йпхярюккхвеяйнцн онкъ  
    call CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD(N,L1,M1,L2,M2,COOAAA,RcoffACF)
    ! опнбепъел бнгмхйюер кх йнлокейярмнярэ йнщттхжхемрю йпхярюкхвеяйнцн онкъ
	DO ISCMME=1,N
       ISA=0
       DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
          IF(RcoffACF(ISCMME,ISA,2).NE.0.D0) THEN
            WRITE(*,*) 'ATTENTION FACTOR OF THE FIELD COMPLEX'
			READ(*,*)
			STOP  
		  ENDIF
       ENDDO
	ENDDO
	
	! тнплхпсел люяяхб пюдхюкэмшу хмрецпюкнб 
	DO ISCMME=1,N
       ISA=0
       DO IKK=IABS(L1-L2),L1+L2,2
          ISA=ISA+1
          RintK(ISCMME,ISA)=CMRI_CENT_R(IKK,COORRR(ISCMME,1),COORRR(ISCMME,2),RO0,H,Npoint,R,Ffun1,Ffun2,RO1)
	   ENDDO
    ENDDO
    
	! бшвхякъел гмювемхе щмепцхх
	EZSUM=0.D0
  	DO ISCMME=1,N
	   ESUM=0.D0
       DO IKK=1,IKNumbre
          ESUM=ESUM+RintK(ISCMME,IKK)*RcoffACF(ISCMME,IKK,1)       
       ENDDO
       EZSUM=EZSUM+Z(ISCMME)*ESUM
	ENDDO
    ! пегскэрюр 
    ENERGYCENTZA=-EZSUM 

    ! сдюкемхе люяяхбнб хг оълърх 
    deallocate(RintK,stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CMME_ONE_ELECTRONIC_INTEGRAL'
      write(*,*) 'THE FILE "RintK" IS NOT REMOVED FROM MEMORY'
	  stop 
    endif
    deallocate(RcoffACF,stat=ierr)
    if(ierr/=0) then
      write(*,*) 'CMME_ONE_ELECTRONIC_INTEGRAL'
      write(*,*) 'THE FILE "RcoffACF" IS NOT REMOVED FROM MEMORY'
	  stop 
    endif
 ENDIF

 !WRITE(8,*) 'ENEW',ENABLA,EKIN,ECENT,ENERGYCENTZA   
 ! пегскэрюр пюяверю люрпхвмнцн щкелемрю
 CMME_ONE_ELECTRONIC_INTEGRAL=ENERGYDIAGONAL+ENERGYCENTZA   
 ! гюохяшбюел йхмерхвеяйсч вюярэ хмрецпюкю
 EnergyKIN=ENABLA+EKIN
 ! гюохяшбюел онремжхюкэмсч вюярэ хмрецпюкю
 EnergyPOT=ECENT+ENERGYCENTZA   


return
end function CMME_ONE_ELECTRONIC_INTEGRAL



! ондопнцпюллю пюяверю люрпхвмнцн щкелемрю дбсущкейрпнммни щмепцхх (йскнмнбяйнцн хмрецпюкю) 
! L1-нпахрюкэмши лнлемр оепбни тсмйжхх
! ML1-опнейжхъ нпахрюкэмнцн лнлемрю оепбни тсмйжхх
! MS1-опнейжхъ яохмю оепбни тсмйжхх
! L2-нпахрюкэмши лнлемр брнпни тсмйжхх
! ML2-опнейжхъ нпахрюкэмнцн лнлемрю брнпни тсмйжхх
! MS2-опнейжхъ яохмю брнпни тсмйжхх
! L3-нпахрюкэмши лнлемр рперэеи тсмйжхх
! ML3-опнейжхъ нпахрюкэмнцн лнлемрю рперэеи тсмйжхх 
! MS3-опнейжхъ яохмю рперэеи тсмйжхх
! L4-нпахрюкэмши лнлемр вербепрни тсмйжхх
! ML4-опнейжхъ нпахрюкэмнцн лнлемрю вербепрни тсмйжхх
! MS4-опнейжхъ яохмю вербепрни тсмйжхх
! H-ьюц
! Npoint-вхякн рнвей
! R(Npoint)-люяяхб гмювемхи юпцслемрю
! Ffun1(Npoint)-люяяхб гмювемхи оепбни тсмйжхх
! Ffun2(Npoint)-люяяхб гмювемхи брнпни тсмйжхх
! Ffun3(Npoint)-люяяхб гмювемхи рперэеи тсмйжхх
! Ffun4(Npoint)-люяяхб гмювемхи вербепрни тсмйжхх
! RO1(Npoint)-люяяхб оепеяверю тсмйжхи (ОЕПБЮЪ ОПНХГБНДМНЪ МНБНИ ОЕПЕЛЕММНИ) 
real(8) function  CMME_TWO_ELECTRONIC_INTEGRAL(L1,ML1,MS1,L2,ML2,MS2,L3,ML3,MS3,L4,ML4,MS4,H,Npoint,R,Ffun1,Ffun2,Ffun3,Ffun4,RO1,FF,F)
 use mcmri,only:CMRI_COULOMB_INTEGRAL
 use mc3js,only:C3JS_Ckq
 implicit none
 integer::L1,ML1,MS1,L2,ML2,MS2,L3,ML3,MS3,L4,ML4,MS4,Npoint
 real(8)::H,FF
 real(8),dimension(:)::R,Ffun1,Ffun2,Ffun3,Ffun4,RO1,F
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer::ierr,IIUR,MIN,MAX
 real(8)::MatrixEleCkq1,MatrixEleCkq2,RCOULOMB,SUMINT
 
 ! гюмскъел оепед пюявернл
 CMME_TWO_ELECTRONIC_INTEGRAL=0.D0

 IF(MS1.EQ.MS3.AND.MS2.EQ.MS4.AND.(ML1+ML2).EQ.(ML3+ML4)) THEN
  

   MIN=MAX0(IABS(L1-L3),IABS(L2-L4))
   MAX=MIN0(IABS(L1+L3),IABS(L2+L4))
   SUMINT=0.D0
   ! жхйк он цюплнмхйюл
   DO IIUR=MIN,MAX,2
       ! нясыеярбкъел пюявер гмювемхи люрпхвмнцн щкелемрю ятепхвеяйни цюплнмхйх
	   MatrixEleCkq1=C3JS_Ckq(L1,ML1,L3,ML3,IIUR,FF,F)
       MatrixEleCkq2=C3JS_Ckq(L4,ML4,L2,ML2,IIUR,FF,F)
	   RCOULOMB=CMRI_COULOMB_INTEGRAL(IIUR,H,Npoint,R,Ffun1,Ffun2,Ffun3,Ffun4,RO1)
	   SUMINT=SUMINT+MatrixEleCkq1*MatrixEleCkq2*RCOULOMB
   ENDDO

   ! пегскэрюр пюяверю
   CMME_TWO_ELECTRONIC_INTEGRAL=SUMINT
 ENDIF
  
return
end function CMME_TWO_ELECTRONIC_INTEGRAL
	
end module mcmme
