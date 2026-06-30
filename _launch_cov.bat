@echo off
cd /d C:\Users\shaya\Documents\VSCode\FSEM\FSEM
set RS="C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
set NW=%1
if "%NW%"=="" set NW=12
for /L %%w in (1,1,%NW%) do start "covw%%w" /min %RS% _cov_worker.R %%w %NW%
echo launched %NW% detached workers
