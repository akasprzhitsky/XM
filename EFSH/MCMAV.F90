      ! THE MODULE OF CALCULATION OF MOLECULAR ANGULAR VALUES VER 1.0 11.2005
	! бяе опюбю гюыхыемш. йнлоюмхъ LabComputerSistemTecnologyKasprzhitsky
	! VER 1.0 NEW  11,2005 цнд

	! лндскэ пюяверю лнкейскъпмшу сцкнбшу гмювемхи бепяхъ 1.0    
	

module mcmav
  implicit none
	
      

	contains



    ! ондопнцпюллю пюяверю сцкнбни ярпсйрспш ноепюрнпю бгюхлндеиярбхъ 
	! ъдеп лнкейскш х щкейрпнмнб (щкейрпнм б йпхярюккхвеяйнл онке)
	! нохяюмхе оюпюлерпнб ондопнцпюллш
   	! N-вхякн ъдеп
	! L1-нпахрюкэмши лнлемр оепбнцн щкейрпнмю
	! ML1-опнейжхъ нпахрюкэмнцн лнлемрю оепбнцн щкейрпнмю
	! L2-нпахрюкэмши лнлемр брнпнцн щкейрпнмю
	! ML2-опнейжхъ нпахрюкэмнцн лнлемрю брнпнцн щкейрпнмю    
	! AC(N,2)-люяяхб сцкнбшу йннпдхмюр ъдеп (б пюдхюмюу) 
	! AC(N,1)-сцнк рерю (0<X<=pi),
	! AC(N,2)-сцнк тх (0<Y<=2*pi) 
	! RcoffCF(N,Numbre,2)-люяяхб сцкнбшу йнщттхжхемрнб йпхярюккхвеяйнцн онкъ 
	! RcoffCF(N,Numbre,1)-пеюкэмюъ вюярэ сцкнбнцн йнщттхжхемрю (онкебнцн) 
	! RcoffCF(N,Numbre,2)-лмхлюъ вюярэ сцкнбнцн йнщттхжхемрю   (онкебнцн) 
	! онъямемхе Numbre-вхякн цюплнмхй
	subroutine CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD(N,L1,ML1,L2,ML2,AC,RcoffCF)
	  use mc3js,only:C3JS_VAR,C3JS_Ckq,C3JS_CkqGG,C3JS_CkqISH
      implicit none
      
	  integer::N,L1,ML1,L2,ML2
      real(8),dimension(:,:)::AC
      real(8),dimension(:,:,:)::RcoffCF
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer,parameter::NmassivVar=800
	  integer ::Numbre,IX,IY,ierr
	  real(8):: FF,MatrixEleCkq,RAS
	  real(8),parameter::PInumbre=3.14159265358979D0
      real(8),allocatable,dimension(:)::Ykq,YSkq,F

      !бшдекъел оюлърэ дкъ люяяхбнб
	  allocate(Ykq(2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "Ykq" IS NOT SELECTED'
	     stop 
	  endif
      allocate(YSkq(2),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "Ykq" IS NOT SELECTED'
	     stop 
	  endif
      allocate(F(NmassivVar),stat=ierr)
	  if(ierr/=0) then
         write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
	     write(*,*) 'MEMORY ON THE FILE "F" IS NOT SELECTED'
	     stop 
	  endif


      ! гюмскъел оепед пюявернл
      RcoffCF=0.D0
	  Ykq=0.D0
	  YSkq=0.D0
      F=0.D0
      ! бяонлнцюрекэмюъ опнжедспю дкъ пюяверю люрпхвмшу щкелемрнб ятепхвеяйни цюплнмхйх
	  !call C3JS_VAR(NmassivVar,FF,F)


      Numbre=0
      ! жхйк он цюплнмхй
      DO IX=IABS(L1-L2),L1+L2,2
	     Numbre=Numbre+1 
	     ! нясыеярбкъел пюявер гмювемхи люрпхвмнцн щкелемрю ятепхвеяйни цюплнмхйх
	     MatrixEleCkq=C3JS_CkqGG(L1,ML1,L2,ML2,IX) !C3JS_Ckq(L1,ML1,L2,ML2,IX,FF,F)
		 !write(6,*) 'ckq',MatrixEleCkq !,C3JS_CkqISH(L1,ML1,L2,ML2,IX) !,C3JS_Ckq(L1,ML1,L2,ML2,IX,FF,F)
	   	 ! йнщттхжхемр
	     RAS=DSQRT(4.d0*PInumbre/float(2*IX+1))
         ! жхйк он мнлепюл ъдеп
	     DO IY=1,N
            ! нясыеярбкъел пюявер гмювемхи ятепхвеяйни цюплнмхйх
            call CMAV_VALUE_SPHERICAL_FUNCTION(IX,ML1-ML2,AC(IY,1),AC(IY,2),Ykq,YSkq)
	        ! пеюкэмюъ вюярэ йнщттхжхемрю
	        RcoffCF(IY,Numbre,1)=MatrixEleCkq*RAS*YSkq(1)
	        ! лмхлюъ вюярэ йющттхжхемрю
            RcoffCF(IY,Numbre,2)=MatrixEleCkq*RAS*YSkq(2)
         ENDDO
      ENDDO


      ! сдюкемхе люяяхбнб хг оълърх   
      deallocate(Ykq,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "Ykq" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(YSkq,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "YSkq" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif
	  deallocate(F,stat=ierr)
	  if(ierr/=0) then
	     write(*,*) 'CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD'
         write(*,*) 'THE FILE "F" IS NOT REMOVED FROM MEMORY'
	     stop 
	  endif


      return
    end subroutine CMAV_ANGULAR_STRUCTURE_CRYSTAL_FIELD










	! ондопнцпюллю пюяверю вхякнбнцн гмювемхъ ятепхвеяйни цюплнмхйх дкъ K х бяеу опнейжхи
    ! нохяюмхе оюпюлерпнб ондопнцпюллш
	! K-мнлеп ятепхвеяйни цюплнмйх
    ! RTeta-СЦНК Б ЯТЕПХВЕЯЙНИ ЯХЯРЕЛШ ЙННПДХМЮР РЕРРЮ (0<RTeta<=pi)
	! RFu-СЦНК Б ЯТЕПХВЕЯЙНИ ЯХЯРЕЛЕ ЙННПДХМЮР ТХ (0<RFu<=2*pi)
	! Ykq(2*K+1,2)-люяяхб гмювемхи ятепхвеяйни цюплнмхйх бяеу опнейжхи дкъ дюммнцн K
	! Ykq(2*K+1,1)-пеюкэмюъ вюярэ ятепхвеяйни цюплнмхйх
	! Ykq(2*K+1,2)-лмхлюъ вюярэ ятепхвеяйни цюплнмхйх 
	! YSkq(2*K+1,2)-люяяхб гмювемхи ятепхвеяйни цюплнмхйх йнлокейямн-янопъфеммни бяеу опнейжхи дкъ дюммнцн K 
	subroutine CMAV_SPHERICAL_FUNCTION(K,RTeta,RFu,Ykq,YSkq)
	  implicit none
      
	  integer::K
	  real(8)::RTeta,RFu
      real(8),dimension(:,:)::Ykq,YSkq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  real(8),parameter::PInum=3.14159265358979D0
	  integer::IZX,IZX1q
      real(8)::RTETAF
      
      ! гюмскъел оепед пюявернл
      Ykq=0.D0
      YSkq=0.D0

      IZX1q=-K
	  ! жхйк он онкнфхрекэмшл опнейжхъл цюплнмхйх бйкчвюъ мнкэ
      do IZX=K+1,2*K+1
         ! опнейжхъ цюплнмхйх
	     IZX1q=IZX1q+(IZX-1)
         ! бшвхякъел гмювемхъ рерю тсмйжхх
         RTETAF=CMAV_VALUE_TETA_FUNCTION_POINT(K,IZX1q,RTeta)
         ! гюохяшбюел цюплнмхйс я мнлепнл K х опнейжхеи IZX1q
	     ! пеюкэмюъ вюярэ 
	     Ykq(IZX,1)=DCOS(float(IZX1q)*RFu)*RTETAF/DSQRT(2.D0*PInum)
         ! лмхлюъ вюярэ
	     Ykq(IZX,2)=DSIN(float(IZX1q)*RFu)*RTETAF/DSQRT(2.D0*PInum)
         ! гюохяшбюел цюплнмхйс я мнлепнл й х опнейжхеи  -IZX1q
         ! пеюкэмюъ вюярэ 
	     Ykq(K+1-IZX1q,1)=(-1.D0)**(IZX1q)*Ykq(IZX,1)
         ! лмхлюъ вюярэ
	     Ykq(K+1-IZX1q,2)=(-1.D0)**(IZX1q+1)*Ykq(IZX,2)
	  enddo

      ! онксвюел йнлокейямн янопъфеммше цюплнмхйх
	  IZX1q=-K
	  do IZX=1,2*K+1
         ! опнейжхъ цюплнмхйх
	     IZX1q=IZX1q+(IZX-1)
         ! пеюкэмюъ вюярэ
	     YSkq(IZX,1)=DCOS(float(2*IZX1q)*RFu)*Ykq(IZX,1)
	     YSkq(IZX,1)=YSkq(IZX,1)+DSIN(float(2*IZX1q)*RFu)*Ykq(IZX,2)
	     ! лмхлюъ вюярэ
         YSkq(IZX,2)=DCOS(float(2*IZX1q)*RFu)*Ykq(IZX,2)
         YSkq(IZX,2)=YSkq(IZX,2)-DSIN(float(2*IZX1q)*RFu)*Ykq(IZX,1)
      enddo
	 
	  return
    end subroutine CMAV_SPHERICAL_FUNCTION


	! ондопнцпюллю пюяверю вхякнбнцн гмювемхъ ятепхвеяйни цюплнмхйх L,M
    ! нохяюмхе оюпюлерпнб ондопнцпюллш
	! L-нпахрюкэмши лнлемр ятепхвеяйни цюплнмйх
    ! M-опнейжхъ нпахрюкэмнцн лнлемрю ятепхвеяйни цюплнмхйх     
    ! RTeta-СЦНК Б ЯТЕПХВЕЯЙНИ ЯХЯРЕЛШ ЙННПДХМЮР РЕРРЮ (0<RTeta<=pi)
	! RFu-СЦНК Б ЯТЕПХВЕЯЙНИ ЯХЯРЕЛЕ ЙННПДХМЮР ТХ (0<RFu<=2*pi)
	! Ykq(2)-люяяхб гмювемхи ятепхвеяйни цюплнмхйх дкъ дюммнцн L х опнейжхх M 
	! Ykq(1)-пеюкэмюъ вюярэ ятепхвеяйни цюплнмхйх
	! Ykq(2)-лмхлюъ вюярэ ятепхвеяйни цюплнмхйх 
	! YSkq(2)-люяяхб гмювемхи ятепхвеяйни цюплнмхйх йнлокейямн-янопъфеммни дкъ дюммнцн K х опнейжхх M 
    subroutine CMAV_VALUE_SPHERICAL_FUNCTION(L,M,RTeta,RFu,Ykq,YSkq)
	  implicit none
      
	  integer::L,M
	  real(8)::RTeta,RFu
      real(8),dimension(:)::Ykq,YSkq
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	  real(8),parameter::PInum=3.14159265358979D0
      real(8)::RTETAF,REY,RIMY
       
      ! гюмскъел оепед пюявернл
      Ykq=0.D0
      YSkq=0.D0

      ! бшвхякъел гмювемхъ рерю тсмйжхх
      RTETAF=CMAV_VALUE_TETA_FUNCTION_POINT(L,IABS(M),RTeta)
      REY=DCOS(float(IABS(M))*RFu)*RTETAF/DSQRT(2.D0*PInum)
      RIMY=DSIN(float(IABS(M))*RFu)*RTETAF/DSQRT(2.D0*PInum)
      ! сярюмюбкхбюел гмюй опнейжхх  
      IF(M.GE.0) THEN
         ! пеюкэмюъ вюярэ 
	     Ykq(1)=REY
         ! лмхлюъ вюярэ
	     Ykq(2)=RIMY
        ELSE
         ! пеюкэмюъ вюярэ 
	     Ykq(1)=(-1.D0)**IABS(M)*REY
         ! лмхлюъ вюярэ
	     Ykq(2)=(-1.D0)**(IABS(M)+1)*RIMY
      ENDIF
      

      ! онксвюел йнлокейямн янопъфеммше цюплнмхйх
      ! пеюкэмюъ вюярэ
      YSkq(1)=DCOS(float(2*M)*RFu)*Ykq(1)
      YSkq(1)=YSkq(1)+DSIN(float(2*M)*RFu)*Ykq(2)
	  ! лмхлюъ вюярэ
      YSkq(2)=DCOS(float(2*M)*RFu)*Ykq(2)
      YSkq(2)=YSkq(2)-DSIN(float(2*M)*RFu)*Ykq(1)
      
      return
    end subroutine CMAV_VALUE_SPHERICAL_FUNCTION


 
   ! ондопнцпюллю онксвемхъ гмювемхъ рерю тсмйжхх (пеййспемрмшл яонянанл)
   ! нохяюмхе оюпюлерпнб ондопнцпюллш
   ! L-нпахрюкэмши лнлемр тсмйжхх
   ! M-опнейжхъ нпахрюкэмнцн лнлемрю (опнцпюллю пюяялюрпхбюер яксвюи (M>=0)  
   ! X-гмювемхе юпцслемрю (б хмрепбюке нр (0,PI))
   real(8) function CMAV_VALUE_TETA_FUNCTION_POINT(L,M,X)
     implicit none
      
	 integer::L,M
     real(8)::X
     !!!!!!!!!!!!!!!!!!!!!!!!!
     integer::IAA2
     real(8)::RHR,Rcoff1,Rcoff2,Rcoff3,RSDXF,RSDXX

     
     ! якюцюелне Teta M,M
	 RSDXF=DSQRT(CMAV_FACTORIAL(2*IABS(M)))/CMAV_FACTORIAL(IABS(M))
	 RSDXX=(DSIN(X*0.5D0)*DCOS(X*0.5D0))**IABS(M)
     RSDXX=RSDXX*RSDXF/DSQRT(2.D0)
	 Rcoff1=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+1))*RSDXX
     !  якюцюелне Teta M+1,M 
	 RSDXF=DSQRT(CMAV_FACTORIAL(2*IABS(M)+1))/CMAV_FACTORIAL(IABS(M))
	 RSDXX=(DSIN(X*0.5D0)*DCOS(X*0.5D0))**IABS(M)
     RSDXX=RSDXX*RSDXF*DCOS(X)/DSQRT(2.D0)
	 Rcoff2=(-1.d0)**IABS(M)*SQRT(float(2*IABS(M)+3))*RSDXX
	 IF(IABS(M).EQ.L) Rcoff3=Rcoff1
     IF(IABS(M).EQ.L-1) Rcoff3=Rcoff2
     ! жхйк он лнлемрюл
	 DO IAA2=IABS(M),L-2
        RSDXF=SQRT(float(2*IAA2+3)*float(2*IAA2+5))
        RSDXF=RSDXF/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
          
	    RSDXX=SQRT(float(2*IAA2+5))
        RSDXX=RSDXX*SQRT(float(IAA2+1-IABS(M))*float(IAA2+1+IABS(M)))
        RSDXX=RSDXX/SQRT(float(2*IAA2+1))
        RSDXX=RSDXX/SQRT(float(IAA2+2-IABS(M))*float(IAA2+2+IABS(M)))
        ! пюявер якедсчыецн щкелемрю
	    Rcoff3=RSDXF*DCOS(X)*Rcoff2-RSDXX*Rcoff1
	    ! ондцнрнбйю й пюяверс якедсчыецн щкелемрю
        Rcoff1=Rcoff2
	    Rcoff2=Rcoff3
     ENDDO
      
	 ! гюохяшбюел гмювемхе тсмйжхх
	 CMAV_VALUE_TETA_FUNCTION_POINT=Rcoff3 

     return
    end function CMAV_VALUE_TETA_FUNCTION_POINT


     
    ! ондопнцпюллю пюяверю тюйрнпхюкю
	! нохяюмхе оюпюлерпнб
    ! NOPI-жекне онкнфхрекэмне вхякн 
	real(8) function CMAV_FACTORIAL(NOPI)
      implicit none

	  integer:: NOPI,IBKL
      real(8):: PROSS
      
	  IF(NOPI.LT.0) THEN
         WRITE(*,*) 'ERROR FACTORIAL (N<0)'
         READ(*,*)
	     STOP
 	  ENDIF

      IF(NOPI.EQ.0) THEN
         CMAV_FACTORIAL=1.D0
         return
	  ENDIF
      
	  PROSS=1.D0
	  do IBKL=1,NOPI
         PROSS=PROSS*float(IBKL)
	  enddo

	  CMAV_FACTORIAL=PROSS

   
      return
    end function CMAV_FACTORIAL
	
end module mcmav
