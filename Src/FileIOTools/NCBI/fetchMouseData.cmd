SETLOCAL ENABLEEXTENSIONS
SET me=%~n0
SET parent=%~dp0
SET fileloc=%parent%SRR_Acc_List_Mouse.txt
SET sradir=C:\Users\dlee\Desktop\SRA\sratoolkit.2.8.2-1-win64\bin
cd %sradir%
for /f "delims=" %%x in (%fileloc%) do (SET srafile=%%x && call :subroutine)
:eof 
PAUSE

:subroutine
prefetch %srafile%
fastq-dump %srafile% --outdir %parent%
GOTO :eof
PAUSE