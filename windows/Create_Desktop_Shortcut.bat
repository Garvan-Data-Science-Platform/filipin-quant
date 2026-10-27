@echo off
REM Double-click this to put a "Filipin Quant" icon on your Desktop.
REM It just runs the PowerShell script with the execution policy bypassed.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0Create_Desktop_Shortcut.ps1"
