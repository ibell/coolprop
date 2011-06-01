mkdir mexBuild
mex -c Nitrogen.c -outdir mexBuild
mex -c Argon.c -outdir mexBuild
mex -c Water.c -outdir mexBuild
mex -c R744.c -outdir mexBuild
mex -c R717.c -outdir mexBuild
mex -c R404A.c -outdir mexBuild
mex -c R410A.c -outdir mexBuild
mex -c R407C.c -outdir mexBuild
mex -c R507A.c -outdir mexBuild
mex -c R290.c -outdir mexBuild
mex -c R134a.c -outdir mexBuild
mex -c R32.c -outdir mexBuild
mex -c Brine.c -outdir mexBuild
mex -c CoolProp.c -outdir mexBuild
mex  mexCoolProp.c -output CoolProp mexBuild/Nitrogen.obj mexBuild/Argon.obj mexBuild/Water.obj mexBuild/R744.obj mexBuild/R717.obj mexBuild/R404A.obj mexBuild/R410A.obj mexBuild/R407C.obj mexBuild/R507A.obj mexBuild/R290.obj mexBuild/R134a.obj mexBuild/R32.obj mexBuild/Brine.obj mexBuild/CoolProp.obj
rmdir('mexBuild','s')