
To Use
======

Get the dll folder from Mathematica by doing::

    FileNameJoin[{$BaseDirectory, "SystemFiles", "LibraryResources", $SystemID}]
    
If this path doesn't exist, make it.

Put the DLL in this folder. (or any other folder in $LibraryPath)

To Build
========

Get "WolframLibrary.h" from C:\\Program Files\\Wolfram Research\\Mathematica\\9.0\\SystemFiles\\IncludeFiles\\C (or similary for your computer) and copy it to the source folder wrappers/Mathematica

