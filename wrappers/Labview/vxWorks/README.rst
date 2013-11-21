Steps
-----

1. Download a WinXP 32-bit version virtual box image
2. Install Winxp 32-bit virtual box image, install 7zip, tortoisegit, git
3. Take the updated_vxworks63gccdist.ZIP from (http://www.ni.com/white-paper/5694/en/) and extract to c:\gccdist so that you have folders c:\gccdist\docs, c:\gccdist\supplemental, c:\gccdist\WindRiver
4. Check out coolprop sources
5. Open a console and cd to CoolProp sources
6. cd to wrappers/Labview/vxWorks
7. run build.bat (It can take a really long time for some reason for Mixtures.cpp, be patient)
8. Upload the generated file CoolProp.out to ni-rt/system folder using the National Instruments Measurement & Automation Explorer.  Right click on the unit, then file transfer. Also see http://www.ni.com/white-paper/3365/en/