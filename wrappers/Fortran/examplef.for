! Example fortran file 

      program EXAMPLE
      double precision T,Q,h,s
      character*20 Ref
      character*2 Output, Name1, Name2
      double precision outVal,Prop1,Prop2

      T = 285
      Q = 0

      Output = 'P'
      Name1  = 'T'
      Prop1  = T
      Name2  = 'Q'
      Prop2  = Q
      Ref    = "R134a"

      Output(LEN(Output):LEN(Output)) = CHAR(0)
      Name1(LEN(Name1):LEN(Name1)) = CHAR(0)
      Name2(LEN(Name2):LEN(Name2)) = CHAR(0)
      Ref(LEN(Ref):LEN(Ref)) = CHAR(0)

      write(*,*) "Saturation pressure for R134a: "
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "Result was: ", outVal/1e2, " bar"

      end program 
