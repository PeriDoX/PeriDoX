::-----------------------------------------------------------------------------
::
:: Purpose:         Create all glossaries for a tex document
::
:: Conventions:
::
::   - Acronyms:    - Input file extension:	    acn*
::                  - Output file extension:    acr*
::                  - Log file extension:       alg*
::
::   - Glossary:    - Input file extension:	    glo*
::                  - Output file extension:    gls*
::                  - Log file extension:       glg*
::
::   - Symbols :    - Input file extension:	    syg*
::                  - Output file extension:    syi*
::                  - Log file extension:       slg*
::
:: Remarks:         - This script has to be physically present in the tex directory.
::                    A symbolic link does not work.
::
:: Author:          Martin Raedel, DLR-FA-STM, WiMi, 10.02.2016
::-----------------------------------------------------------------------------
::
:: @echo off verhindert Ausgabe Befehlszeilen bis Stapelverarbeitung beendet
@echo off
::
::-----------------------------------------------------------------------------
:: set variables - no spaces between variable name and = allowed
::-----------------------------------------------------------------------------
:: 
:: manual mode:
::set maintexname=Mixed-Dimensional-Coupling
::
:: automatic mode:
set glo_num=2
::
set glo_input_list[0]=acn
set glo_input_list[1]=glo
set glo_input_list[2]=syg
::
set glo_output_list[0]=acr
set glo_output_list[1]=gls
set glo_output_list[2]=syi
::
set glo_log_list[0]=alg
set glo_log_list[1]=glg
set glo_log_list[2]=slg
::
::-----------------------------------------------------------------------------
:: create glossaries
::-----------------------------------------------------------------------------
::
:: Acronyms:
::
::makeindex -s %maintexname%.ist -t %maintexname%.alg -o %maintexname%.acr %maintexname%.acn
::
:: Glossar:
::
::makeindex -s %maintexname%.ist -t %maintexname%.glg -o %maintexname%.gls %maintexname%.glo
::
:: Symbols:
::   manual:
::makeindex -s %maintexname%.ist -t %maintexname%.slg1 -o %maintexname%.syi1 %maintexname%.syg1
::   with automatic symbol input file search
setlocal EnableDelayedExpansion
::
FOR /L %%i IN (0, 1, %glo_num%) DO (
  set inputextension=!glo_input_list[%%i]!
  set outputextension=!glo_output_list[%%i]!
  set logextension=!glo_log_list[%%i]!
  ::echo "inputextension: !inputextension!"
  FOR %%f IN (*.!inputextension!*) DO (
    ::echo %%f
    set fullfilename=%%f
    set filename=%%~nf
    set fileextensionwithdot=%%~xf
    set fileextension=!fileextensionwithdot:~1!
    ::echo "fullfilename: !fullfilename! | fileextension: !fileextension! | !fileextension:~3,1!"
    IF DEFINED fileextension IF "!fileextension:~3,1!"=="" (
      ::echo 3 or less characters
      makeindex -s !filename!.ist -t !filename!.!logextension! -o !filename!.!outputextension! !filename!!fileextensionwithdot! >nul 2>&1
    ) ELSE (
      ::echo more than 3 characters
      set index=!fileextension:~3,1!
      ::echo !index!
      makeindex -s !filename!.ist -t !filename!.!logextension!!index! -o !filename!.!outputextension!!index! !filename!!fileextensionwithdot!>nul 2>&1
    )
  )
)
::
endlocal