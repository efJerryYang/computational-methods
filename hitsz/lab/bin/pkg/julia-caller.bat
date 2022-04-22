@echo off

title Computational Method Lab Demo
color f
set juliapath=%~dp0bin\
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
color f
chcp 65001
echo:
echo        Welcome! 
echo:       
echo        This is Computational Method Lab Demo, which is used to demonstrate the results of the 5
echo        lab questions.
echo:       
goto menu

:exit_info
echo:
echo Press any key to continue... & pause > nul
goto end

:menu
color f
echo:
echo        Please enter the number according to the prompt, run the corresponding lab code or exit:
@REM echo        请根据提示输入数字, 运行相应实验代码或退出：
echo:
echo         0. Run All Labs
echo         1. Lagrange
echo         2. Romberg
echo         3. Runge-Kutta
echo         4. Newton
echo         5. Gauss
echo        -1. Exit
echo:
set num=-999
set /p num=LAB-DEMO# 
if %num%==0 ( 
    echo:
    echo Running all the labs...
    goto runall
) else if %num%==1 (
    echo:
    echo Lab %num%: Lagrange
    goto lagrange
) else if %num%==2 (
    echo:
    echo Lab %num%: Romberg
    goto romberg
) else if %num%==3 ( 
    echo:
    echo Lab %num%: Runge-Kutta
    goto rungekutta
) else if %num%==4 (
    echo:
    echo Lab %num%: Newton
    goto newton
) else if %num%==5 ( 
    echo:
    echo Lab %num%: Gauss
    goto gauss
) else if %num%==-1 (
    goto exit_info
) else (
    echo:
    echo Please enter a valid number and try again.
    echo:
    @REM chcp 65001
    goto greeting
)

@REM :start_fail
@REM echo:
@REM echo 错误! 无法启动程序 julia.exe.
@REM @REM echo ERROR! Fail to start julia.exe.
@REM echo:
@REM if exist "%~dp0pkg" ( 
@REM     echo 文件夹 pkg 存在, 为便于重新安装请将其所有子目录移除.
@REM     @REM echo Folder pkg exists, please remove it and all its subdirectories.
@REM ) else (
@REM     echo 文件夹 pkg 不存在.
@REM     @REM echo Folder pkg not exists.
@REM )
@REM echo:
@REM echo 请将 pkg.zip 解压到 .\pkg 目录下, 并确保文件夹 pkg 没有出现冗余的嵌套.
@REM echo:
@REM echo 当前的工作目录树形结构如下:
@REM echo:
@REM type .\etc\folder.txt
@REM echo:
@REM echo:
@REM echo 如果有需要, 请在 .\etc\file.txt 中获取更多有关依赖的信息.
@REM goto exit_info

:runall
echo:
echo    It really takes some time to compile the program at runtime...
echo:
call "%juliapath%Demo.exe" -a
echo:
goto greeting

:lagrange
echo:
echo    It really takes some time to compile the program at runtime...
echo:
call "%juliapath%Demo.exe" -1
@REM call "%juliapath%julia.exe" --sysimage=sys_lab.dll -q -- lab1-lagrange.jl
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
@REM chcp 65001
goto greeting

:romberg
echo:
echo    It really takes some time to compile the program at runtime...
echo:
call "%juliapath%Demo.exe" -2
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
@REM chcp 65001
goto greeting

:rungekutta
echo:
echo    It really takes some time to compile the program at runtime...
echo:
call "%juliapath%Demo.exe" -3
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
@REM chcp 65001
goto greeting

:newton
echo:
echo    It really takes some time to compile the program at runtime...
echo:
call "%juliapath%Demo.exe" -4
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
@REM chcp 65001
goto greeting

:gauss
call "%juliapath%Demo.exe" -5
echo:
echo    It really takes some time to compile the program at runtime...
echo:
@REM color f
echo:
@REM echo Press any key to continue... & pause > nul
@REM chcp 65001
goto greeting

:end
exit