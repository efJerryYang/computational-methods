@echo off
rem for and del

for /d %%a in (julia-1.7.2/*) do if  %%a==bin (rd julia-1.7.2/%%a)

pause > nul