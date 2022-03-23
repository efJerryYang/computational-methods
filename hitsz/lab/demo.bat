@echo off
chcp 65001
title Lab Demo: Computational Method
echo:
echo Press any key to continue...
pause > nul

set juliapath=%~dp0pkg\bin\julia.exe

if exist %juliapath% ( goto menu ) else ( goto start_fail )

:exit_info
echo:
echo If some characters cannot be displayed correctly, please make sure the terminal encoding is UTF-8.
echo:
echo That means you should see "Active code page: 65001" at the entry of this terminal.
echo:
echo Press any key to continue...
pause > nul
goto end

:menu
echo:
echo 请根据提示输入数字, 运行相应实验代码或退出
echo:
echo   1. Lagrange
echo   2. Romberg
echo   3. Runge-Kutta
echo   4. Newton
echo   5. Gauss
echo   0. Exit
echo:
set num=-1
set /p num=请输入数字:
if %num%==1 ( 
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
    goto runge_kutta
) else if %num%==4 (
    echo:
    echo Lab %num%: Newton
    goto newton
) else if %num%==5 ( 
    echo:
    echo Lab %num%: Gauss
    goto gauss
) else if %num%==0 (
    goto exit_info
) else (
    echo:
    echo 请输入有效的数字重试.
    goto menu
)

@REM echo 解释代码需要一些时间...
@REM echo It really takes time to run the program...
@REM echo:
@REM call %juliapath% %~dp0labs\lab2.jl
@REM call %juliapath%
goto exit_info


:start_fail
echo 错误! 无法启动程序 julia.exe.
@REM echo ERROR! Fail to start julia.exe.
echo:
if exist %~dp0pkg ( 
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
goto menu

:romberg
goto menu

:runge_kutta
goto menu

:newton
goto menu

:gauss
goto menu

:end