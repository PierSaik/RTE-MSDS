@echo off
cls

echo #=======================================================#
echo #                                                       #              
echo #     RECHERCHE DU DOMAINE DE SIMULATION STABLE         #
echo #                                                       #
echo #=======================================================#
echo.
echo #=======================INFORMATION=====================#
echo #                                                       #
echo # Vous devez utiliser ce programme dans le repertoire   #
echo # contenant un fichier .dta et les regulateurs associes #
echo # aux machine ainsi qu'un dictionnaire                  #
echo #                                                       #
echo #=======================================================#
echo.

set /p RefdtaFile= Indiquer le fichier .dta : 
echo.
set /p generatorsFile= Indiquer le fichier contenant les generateurs a traiter : 
echo.
set /p dictFile= Indiquer le fichier dictionnaire : 
echo.

echo Les generateurs a traite sont :
type %generatorsFile%
echo.

echo.
echo Confirmer ?
pause

python ./IndusMachineDomain.py %RefdtaFile% %generatorsFile% %dictFile%

pause