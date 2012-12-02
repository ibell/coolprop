Ian Bell, Ph.D. (ian.h.bell@gmail.com)
December 2, 2012

Info
----
These are the files needed for Csharp applications on Windows. A similar build process should be used for 
non-Windows application though no knowledge is available for non-Windows applications

Build
-----
The build batch file BuildCsharpDLL.bat should be run to generate the C# files (they are put into the VSCsharp folder)

The Visual Studio C# file VSCsharp/CoolProp.NET.csproj can be run to generate the executable.  Make sure you keep it in "x86" 
rather than "Any CPU" architecture so that it uses 32-bit calling conventions.  Otherwise you will get PINVOKE errors!!!
