@echo off
REM Deploy the "Filipin Quant" icon to the shared Public Desktop, so it appears
REM for EVERY user of this machine (e.g. installing for another user from an admin account).
REM Right-click -> "Run as administrator".
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0Create_Desktop_Shortcut.ps1" -AllUsers
