@echo off

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

rem ----------------------------------------------------------------------
rem ifort compile
rem ----------------------------------------------------------------------
ifort .\src\iric.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort  /F100000000 /f77rtl /fixed /extend-source:132 /MD /c .\src\nays2d+.f
ifort *.obj .\lib\iriclib.lib -o ".\install\nays2d+.exe"

del *.obj
del *.mod
