@echo off
title lab 1
@REM pause

@REM @REM :label1
@REM @REM echo "Hello world! This is label 1."
@REM @REM goto label3

@REM :label2
@REM echo Hello world! This is label 2.
@REM pause

@REM set clk=10
@REM echo %clk%
@REM pause
@REM @REM echo key=
@REM set /p key=input key number:
@REM echo key:  %key%
@REM if "%key%"="10" echo good number!
@REM else echo bad number...
@REM pause
@REM exit

@REM @REM :label3
@REM @REM echo This is label 3.
@REM @REM goto label2

@REM dir *julia.exe* /s/b 
dir /s/b 
set v=hello
if %v%==hello (echo ok) else (echo no)
if exist path (echo ok) else (echo no)
pause 
@REM REM  nul