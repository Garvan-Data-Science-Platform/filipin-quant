@echo off
REM ============================================================
REM  Filipin Quantification GUI - Windows launcher
REM  Double-click (or use the desktop shortcut) to start the GUI.
REM
REM  Runs the environment's python.exe DIRECTLY, with the env + conda base
REM  folders on PATH. It does NOT call `conda run` / `activate.bat` (those
REM  recurse and abort on machines with a looping conda AutoRun/hook).
REM
REM  The env lives under C:\ProgramData (all-users, READ-ONLY for a normal
REM  user), so we redirect numba/matplotlib caches to a per-user writable
REM  folder - otherwise numba fails with PermissionError writing into the
REM  env's site-packages, and the GUI never opens for a non-admin user.
REM
REM  Point it at the env ONE of two ways:
REM   - ENV_PREFIX : full path to the env folder (wins if set), e.g.
REM                  C:\ProgramData\miniconda3\envs\cell_quant_sam
REM   - ENV_NAME   : env name; looked up as <conda>\envs\<ENV_NAME> under the
REM                  usual install roots (which also gives us the base folder).
REM ============================================================

set "ENV_PREFIX="
set "ENV_NAME=cell_quant_sam"
set "LOG=%TEMP%\filipin_gui.log"

REM ---- project root = the folder ABOVE this script -----------
pushd "%~dp0.."
set "PROJECT_DIR=%CD%"

REM ---- locate the env's python.exe and its conda base (no conda calls) ----
set "PREFIX="
set "CONDA_BASE="
if defined ENV_PREFIX if exist "%ENV_PREFIX%\python.exe" (
    set "PREFIX=%ENV_PREFIX%"
    for %%i in ("%ENV_PREFIX%\..\..") do set "CONDA_BASE=%%~fi"
)
if not defined PREFIX (
    for %%D in (
        "%PROGRAMDATA%\miniconda3"
        "%PROGRAMDATA%\anaconda3"
        "%USERPROFILE%\miniconda3"
        "%USERPROFILE%\anaconda3"
        "%LOCALAPPDATA%\miniconda3"
        "%LOCALAPPDATA%\anaconda3"
        "C:\miniconda3"
        "C:\anaconda3"
    ) do (
        if not defined PREFIX if exist "%%~D\envs\%ENV_NAME%\python.exe" ( set "PREFIX=%%~D\envs\%ENV_NAME%" & set "CONDA_BASE=%%~D" )
    )
)

echo ============================================================
echo   Filipin Quantification GUI
echo   env:  %PREFIX%
echo   base: %CONDA_BASE%
echo   log:  %LOG%
echo ============================================================
echo.

if not defined PREFIX (
    echo [ERROR] Could not find the env's python.exe.
    echo         Set ENV_PREFIX at the top of this file to the env folder.
    echo.
    pause
    popd
    exit /b 1
)

set "PATH=%PREFIX%;%PREFIX%\Library\mingw-w64\bin;%PREFIX%\Library\usr\bin;%PREFIX%\Library\bin;%PREFIX%\Scripts;%PREFIX%\bin;%CONDA_BASE%;%CONDA_BASE%\Library\bin;%CONDA_BASE%\Scripts;%CONDA_BASE%\condabin;%PATH%"

REM ---- per-user writable cache dirs (the env folder is read-only) ----
set "CACHE=%LOCALAPPDATA%\FilipinQuant"
set "NUMBA_CACHE_DIR=%CACHE%\numba"
set "MPLCONFIGDIR=%CACHE%\matplotlib"
if not exist "%NUMBA_CACHE_DIR%" mkdir "%NUMBA_CACHE_DIR%" 2>nul
if not exist "%MPLCONFIGDIR%" mkdir "%MPLCONFIGDIR%" 2>nul

echo Starting... a napari window will open shortly.
echo (First segmentation downloads model weights - needs internet once.)
"%PREFIX%\python.exe" "%PROJECT_DIR%\filipin_gui.py" > "%LOG%" 2>&1
set "RC=%ERRORLEVEL%"

if not "%RC%"=="0" (
    echo.
    echo [ERROR] The GUI exited with code %RC%. Details:
    echo ------------------------------------------------------------
    type "%LOG%"
    echo ------------------------------------------------------------
    pause
)

popd
exit /b %RC%
