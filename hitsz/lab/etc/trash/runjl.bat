@echo off
@REM echo 
if exist "%~dp0\julia-1.7.2\bin" echo "ehll"
for /r "%~dp0\julia-1.7.2\bin" %%v in (*.exe) do (
    echo %%v
    @REM if "%%v==julia.exe" (echo %%v)
    call %%v
)
pause