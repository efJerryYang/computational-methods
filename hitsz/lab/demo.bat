@echo off
if "%1" == "h" goto begin
mshta vbscript:createobject("wscript.shell").run("""%~nx0"" h",0)(window.close)&&exit
:begin
@echo off

    setlocal enableextensions disabledelayedexpansion

    set "consoleName=executor"

    :: http://technet.microsoft.com/en-us/library/cc978570.aspx
    (   reg add "HKCU\Console\%consoleName%" /f 
        reg add "HKCU\Console\%consoleName%" /f /v "FaceName"                /t "REG_SZ"     /d "Consolas"
        reg add "HKCU\Console\%consoleName%" /f /v "FontFamily"              /t "REG_DWORD"  /d 0x00000036
        reg add "HKCU\Console\%consoleName%" /f /v "FontSize"                /t "REG_DWORD"  /d 0x00120007
        reg add "HKCU\Console\%consoleName%" /f /v "FontWeight"              /t "REG_DWORD"  /d 0x00000001
        reg add "HKCU\Console\%consoleName%" /f /v "QuickEdit"               /t "REG_DWORD"  /d 0x00000001
        reg add "HKCU\Console\%consoleName%" /f /v "ScreenBufferSize"        /t "REG_DWORD"  /d 0x07d00050
        reg add "HKCU\Console\%consoleName%" /f /v "HistoryBufferSize"       /t "REG_DWORD"  /d 0x00000999
        reg add "HKCU\Console\%consoleName%" /f /v "NumberOfHistoryBuffers"  /t "REG_DWORD"  /d 0x00000005
        reg add "HKCU\Console\%consoleName%" /f /v "WindowSize"              /t "REG_DWORD"  /d 0x00200040
        reg add "HKCU\Console\%consoleName%" /f /v "FullScreen"              /t "REG_DWORD"  /d 0x00000001
    ) > nul

    start "%consoleName%" /max "%~dp0src\julia-caller.bat"