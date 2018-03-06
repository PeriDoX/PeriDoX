::-----------------------------------------------------------------------------
:: Delete all files except for necessary LaTeX-files in the root and subdirs
::
:: This script has to be physically present in the tex directory. A symbolic
:: does not work
::
:: Martin Raedel, DLR-FA-STM, WiMi, 10.02.2016
::-----------------------------------------------------------------------------
::
:: @echo off verhindert Ausgabe Befehlszeilen bis Stapelverarbeitung beendet
@echo off
::
::-----------------------------------------------------------------------------
:: set variables - no spaces between variable name and = allowed
::-----------------------------------------------------------------------------
REM ::~~~~~~~~~~~~~~~~~
REM :: DELETE ALL EXCEPT-MODE
REM ::~~~~~~~~~~~~~~~~~
REM :: Anzahl Ordner
REM set keep_num=2
REM :: Ordnerliste
REM set keep_list[0]=cmd
REM set keep_list[1]=tex
REM set keep_list[2]=kilepr

REM setlocal EnableDelayedExpansion

REM echo Main directory: %CD%

REM echo   File attributes for files to keep - Modify
REM for /L %%i IN (0, 1, %keep_num%) DO (
  REM dir %CD%\*.!keep_list[%%i]! | find "<SYMLINK>" && (
    REM echo %CD%\*.!keep_list[%%i]! is a link
  REM )
  REM ATTRIB +r +s %CD%\*.!keep_list[%%i]!
REM )

REM echo   Delete files in %CD%
REM DEL *.* /Q

REM echo   Delete backup files %CD%-subdirectories
REM DEL *.backup /S /Q

REM echo   File attributes for files to keep - Reestablish
REM for /L %%i IN (0, 1, %keep_num%) DO (
  REM ATTRIB -r -s *.!keep_list[%%i]!
REM )

REM endlocal

::~~~~~~~~~~~~~~~~~
:: DELETE SPECIFIC - DUE TO SAFETY ISSUES
::~~~~~~~~~~~~~~~~~
:: Anzahl Ordner
set delete_maindir_num=35
set delete_subdir_num=2
set delete_dir_num=1
:: Ordnerlisten
set delete_maindir_list[0]=acn
set delete_maindir_list[1]=acr
set delete_maindir_list[2]=alg
set delete_maindir_list[3]=aux
set delete_maindir_list[4]=bbl
set delete_maindir_list[5]=bcf
set delete_maindir_list[6]=blg
set delete_maindir_list[7]=dpth
set delete_maindir_list[8]=dvi
set delete_maindir_list[9]=glg
set delete_maindir_list[10]=glo
set delete_maindir_list[11]=gls
set delete_maindir_list[12]=idx
set delete_maindir_list[13]=ilg
set delete_maindir_list[14]=ind
set delete_maindir_list[15]=ist
set delete_maindir_list[16]=lock
set delete_maindir_list[17]=lof
set delete_maindir_list[18]=log
set delete_maindir_list[19]=lol
set delete_maindir_list[20]=lot
set delete_maindir_list[21]=md5
set delete_maindir_list[22]=mw
set delete_maindir_list[23]=nlo
set delete_maindir_list[24]=out
::set delete_maindir_list[14]=pdf
set delete_maindir_list[25]=ps
set delete_maindir_list[26]=toc
set delete_maindir_list[27]=run.xml
set delete_maindir_list[28]=slg
set delete_maindir_list[29]=syg
set delete_maindir_list[30]=syi
set delete_maindir_list[31]=synctex*
set delete_maindir_list[32]=synctex.gz
set delete_maindir_list[33]=tex.backup
set delete_maindir_list[34]=user.adi
set delete_maindir_list[35]=wrt

set delete_subdir_list[0]=aux
set delete_subdir_list[1]=backup
set delete_subdir_list[2]=log

set dir_list[0]=ZZZ_TikZ
set dir_list[1]=ZZZ_Table


:: Action
setlocal EnableDelayedExpansion

echo Main directory: %CD%

echo   Delete files in main-directory
for /L %%i IN (0, 1, %delete_maindir_num%) DO (
  IF EXIST %CD%\*.!delete_maindir_list[%%i]! ( DEL %CD%\*.!delete_maindir_list[%%i]! /Q )
)

echo   Delete files in subdirectories
for /L %%i IN (0, 1, %delete_subdir_num%) DO (
  ::IF EXIST %CD%\*.!delete_subdir_list[%%i]! ( DEL %CD%\*.!delete_subdir_list[%%i]! /Q )
  DEL *.!delete_subdir_list[%%i]! /S /Q
)

endlocal

echo Delete temporary folders

setlocal EnableDelayedExpansion
for /L %%i IN (0, 1, %delete_dir_num%) DO (
  echo   Directory: "!dir_list[%%i]!"
  ::IF EXIST !dir_list[%%i]! ( echo     Exists) ELSE ( echo      Does not exist && exit)
  IF EXIST !dir_list[%%i]! ( 
    echo     Exists
    echo     Delete
    RMDIR /S /Q !dir_list[%%i]!
  ) ELSE ( echo     Does not exist)
)
endlocal