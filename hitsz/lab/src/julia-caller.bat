@REM 不要在代码块内部进行CLS，最好让调用者进行
@echo off
@REM cd ..
title Computational Method Lab Demo
@REM color f
set juliapath=%~dp0
set juliapath=%juliapath:~0,-4%pkg\bin\julia.exe
@REM :~0,-4%pkg\julia.bat
@REM echo %juliapath%
@REM pause
:copyright
cls
chcp 65001
echo:
echo        Computational Method Lab Demo  Copyright (C) 2022  Jerry Yang
echo:       
echo        This program comes with ABSOLUTELY NO WARRANTY; 
echo        for details type 'show -w'.
echo:       
echo        This is free software, and you are welcome to redistribute it under certain conditions;
echo        type 'show -c' for details.
echo:       
echo        View the complete TERMS AND CONDITIONS, type 'show -t'
echo:
echo        If you have any questions about the program, please email me at efjerryyang@outlook.com.
echo:
echo        Just for using the program, press [Enter] key to continue...
echo:
set choice=
set /p choice=LAB-DEMO# 
@REM echo %choice%
if "%choice%"=="" ( 
    cls
    goto greeting
    )
if /i "%choice%"=="show -w" ( 
    cls
    goto warranty 
) else if /i "%choice%"=="show -c" ( 
    cls
    goto conditions 
) else if /i "%choice%"=="show -t" ( 
    cls
    goto terms 
) else ( 
    cls
    goto copyright
    )
pause
cls
goto greeting

:warranty
echo:
echo Disclaimer of Warranty:
echo:
echo        THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
echo        EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE 
echo        THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, 
echo        BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
echo        PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. 
echo        SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR 
echo        OR CORRECTION.
echo:
echo Press any key to continue... & pause > nul
cls
goto copyright

:conditions
echo Basic Permissions:
echo:
echo        All rights granted under this License are granted for the term of copyright on the Program, 
echo        and are irrevocable provided the stated conditions are met. This License explicitly affirms 
echo        your unlimited permission to run the unmodified Program. The output from running a covered 
echo        work is covered by this License only if the output, given its content, constitutes a covered 
echo        work. This License acknowledges your rights of fair use or other equivalent, as provided by 
echo        copyright law.
echo:
echo        You may make, run and propagate covered works that you do not convey, without conditions so 
echo        long as your license otherwise remains in force. You may convey covered works to others for 
echo        the sole purpose of having them make modifications exclusively for you, or provide you with 
echo        facilities for running those works, provided that you comply with the terms of this License 
echo        in conveying all material for which you do not control copyright. Those thus making or 
echo        running the covered works for you must do so exclusively on your behalf, under your direction 
echo        and control, on terms that prohibit them from making any copies of your copyrighted material 
echo        outside their relationship with you.
echo:
echo        Conveying under any other circumstances is permitted solely under the conditions stated below. 
echo        Sublicensing is not allowed; section 10 makes it unnecessary.
echo:
echo Press any key to continue... & pause > nul
cls
goto copyright

:terms
echo:
type "%~dp0LICENSE"
echo:
echo Press any key to continue... & pause > nul
cls
goto copyright

:greeting
@REM color f
chcp 65001
echo:
echo        Welcome! 
echo:       
echo        This is Computational Method Lab Demo, which is used to demonstrate the results of 
echo        the 5 lab questions.
echo:       
@REM echo Press any key to continue... & pause > nul
goto menu

:install
echo    Installing Julia Packages...
set JULIA_PKG_SERVER=mirrors.tuna.tsinghua.edu.cn/julia
call "%juliapath%" "%~dp0julia\install-pkg.jl"
call "%juliapath%" "%~dp0julia\build-pkg.jl"
goto menu
@REM @echo on
@REM set juliapath=%~dp0pkg\bin\julia.exe
@REM 
@REM if exist %juliapath% ( goto menu ) else ( goto start_fail )

:exit_info
echo:
@REM echo If some characters cannot be displayed correctly, please make sure the terminal encoding is UTF-8.
@REM echo:
@REM echo That means you should see "Active code page: 65001" at the entry of this terminal.
@REM echo:
echo Press any key to continue... & pause > nul
goto end

@REM :menu
@REM :MoveCursorUp
@REM setlocal
@REM for /F %%a in ('echo prompt $H ^| cmd') do set "BS=%%a"
@REM for /F "skip=4 delims=pR tokens=2" %%a in ('reg query hkcu\environment /v temp' ) do set "TAB=%%a"
@REM for /F "tokens=2 delims=0" %%a in ('shutdown /? ^| findstr /BC:E') do if not defined TAB set "TAB=%%a"
@REM timeout /T 1 > CON | cmd /Q /C for /F %%C in ('copy /Z "%~F0" NUL') do set /P "=.%%C%TAB%%TAB%%TAB%%TAB%%TAB%%TAB%%TAB%%TAB%" ^& set /P "=.%TAB%%BS%%BS%%%C"

@REM :XY
@REM set "X=%~1"
@REM set "Y=%~2"
@REM set spc=

@REM for /l %%A IN (1,1,!X!) Do Set "spc=!spc! "
@REM for /l %%A IN (1,1,!Y!) Do Echo. 

@REM @REM Goto :eof

:menu
@REM color f
echo:
echo        Please enter the number according to the prompt, run the corresponding lab code or exit:
@REM echo        请根据提示输入数字, 运行相应实验代码或退出：
echo:
echo        1. Lagrange
echo        2. Romberg
echo        3. Runge-Kutta
echo        4. Newton
echo        5. Gauss
echo        0. Exit
echo:
set num=-1
set /p num=LAB-DEMO# 
if %num%==1 ( 
    @REM cls
    echo:
    echo Lab %num%: Lagrange
    goto lagrange
) else if %num%==2 (
    @REM cls
    echo:
    echo Lab %num%: Romberg
    goto romberg
) else if %num%==3 ( 
    @REM cls
    echo:
    echo Lab %num%: Runge-Kutta
    goto rungekutta
) else if %num%==4 (
    @REM cls
    echo:
    echo Lab %num%: Newton
    goto newton
) else if %num%==5 ( 
    @REM cls
    echo:
    echo Lab %num%: Gauss
    goto gauss
) else if %num%==0 (
    goto exit_info
) else (
    @REM cls
    echo:
    echo 请输入有效的数字重试.
    echo:
    chcp 65001
    goto menu
)

@REM echo 解释代码需要一些时间...
@REM echo It really takes time to run the program...
@REM echo:
@REM call %juliapath% %~dp0labs\lab2.jl
@REM call %juliapath%
goto exit_info


:start_fail
echo:
echo 错误! 无法启动程序 julia.exe.
@REM echo ERROR! Fail to start julia.exe.
echo:
if exist "%~dp0pkg" ( 
    echo 文件夹 pkg 存在, 为便于重新安装请将其所有子目录移除.
    @REM echo Folder pkg exists, please remove it and all its subdirectories.
) else (
    echo 文件夹 pkg 不存在.
    @REM echo Folder pkg not exists.
)
echo:
echo 请将 pkg.zip 解压到 .\pkg 目录下, 并确保文件夹 pkg 没有出现冗余的嵌套.
echo:
echo 当前的工作目录树形结构如下:
echo:
type .\etc\folder.txt
echo:
echo:
echo 如果有需要, 请在 .\etc\file.txt 中获取更多有关依赖的信息.
goto exit_info

:lagrange
call "%juliapath%" -q --sysimage="%~dp0sys_test_roots.so" "%~dp0julia\lab1-lagrange.jl"
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
chcp 65001
goto menu

:romberg
call "%juliapath%" -q --sysimage="%~dp0sys_test_roots.so" "%~dp0julia\lab2-romberg.jl"
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
chcp 65001
goto menu

:rungekutta
call "%juliapath%" -q --sysimage="%~dp0sys_test_roots.so" "%~dp0julia\lab3-rungekutta.jl"
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
chcp 65001
goto menu

:newton
call "%juliapath%" -q --sysimage="%~dp0sys_test_roots.so" "%~dp0julia\lab4-newton.jl"
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
chcp 65001
goto menu

:gauss
call "%juliapath%" -q --sysimage="%~dp0sys_test_roots.so" "%~dp0julia\lab5-gauss.jl"
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
chcp 65001
goto menu

:end
exit