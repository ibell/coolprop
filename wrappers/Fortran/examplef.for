! Example fortran file 

      program EXAMPLE
      double precision T,Q,D,h,s
      character(LEN=32) Ref
      character(LEN=2) Output, Name1, Name2
      double precision outVal,Prop1,Prop2

      T = 285
      Q = 0
      D = 1250;

      Output = "P"//CHAR(0)
      Name1  = "T"//CHAR(0)
      Prop1  = T
      Name2  = "Q"//CHAR(0)
      Prop2  = Q
      Ref    = "R134a"//CHAR(0)

!       Output(LEN(Output):LEN(Output)) = CHAR(0)
!       Name1(LEN(Name1):LEN(Name1)) = CHAR(0)
!       Name2(LEN(Name2):LEN(Name2)) = CHAR(0)
!       Ref(LEN(Ref):LEN(Ref)) = CHAR(0)

      write(*,*) "Saturation pressure for R134a: "
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "Result was: ", outVal/1e2, " bar"

      Output = "S"//CHAR(0)

      outVal = props1(Ref , "Tcrit"//CHAR(0))
      write(*,*) "double for props1    : ", outVal 
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "double for props     : ", outVal
      outVal = derivterms("dpdrho"//CHAR(0), Prop1, D, Ref)
      write(*,*) "double for deriv     : ", outVal
      outVal = set_reference_states(Ref, "IIR"//CHAR(0))
      write(*,*) "int for reference    : ", outVal
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "double for props     : ", outVal
      outVal = set_reference_states(Ref, "ASHRAE"//CHAR(0))
      write(*,*) "int for reference    : ", outVal
      outVal = props(Output, Name1, Prop1, Name2,Prop2, Ref)
      write(*,*) "double for props     : ", outVal
      outVal = enable_ttse_lut(Ref)
      write(*,*) "bool for enable ttse : ", outVal
      outVal = isenabled_ttse_lut(Ref)
      write(*,*) "bool for check ttse  : ", outVal
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "double for props     : ", outVal
      outVal = set_ttse_mode(Ref , "bicubic"//CHAR(0))
      write(*,*) "int for set ttse mode: ", outVal
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "double for props     : ", outVal
      outVal = disable_ttse_lut(Ref)
      write(*,*) "bool for disable ttse: ", outVal
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "double for props     : ", outVal

      end program 
