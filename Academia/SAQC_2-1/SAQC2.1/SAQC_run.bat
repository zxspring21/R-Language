@echo off
mode 125,40 
set yourDir=C:\Program Files\R\
If Not exist "C:\Program Files\R"\ (
call %1%R-3.0.1-win.exe
if %PROCESSOR_ARCHITECTURE%==x86 (
set List=C:\Program Files\R\R-3.0.1\bin\i386\
) else (
set List=C:\Program Files\R\R-3.0.1\bin\x64\
)
) else (
IF %PROCESSOR_ARCHITECTURE%==x86 (
for /d %%i in ("%yourDir%*") do (
If Not exist %%~i\bin\ (
call %1%R-3.0.1-win.exe
set List=C:\Program Files\R\R-3.0.1\bin\i386\
) else (
set List=%%~i\bin\i386\
)
)) else (
for /d %%i in ("%yourDir%*") do (
If Not exist %%~i\bin\ (
call %1%R-3.0.1-win.exe
set List=C:\Program Files\R\R-3.0.1\bin\x64\
) else (
set List=%%~i\bin\x64\
)
))
)
set Path=%Path%;%List%
R --no-save --slave < SAQC_path.R