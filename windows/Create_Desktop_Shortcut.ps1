# Creates a "Filipin Quant" shortcut that launches the GUI.
#
#   (default)      -> current user's Desktop
#   -AllUsers      -> the shared Public Desktop (visible to every user,
#                     e.g. deploying for someone else; needs admin rights)
#
# Run once. Right-click -> Run with PowerShell, or use the .bat wrappers.

param([switch]$AllUsers)

$here     = Split-Path -Parent $MyInvocation.MyCommand.Definition
$target   = Join-Path $here "launch_filipin_gui.bat"
$icon     = Join-Path $here "filipin.ico"
$workdir  = Split-Path -Parent $here            # project root

if ($AllUsers) {
    $desktop = Join-Path $env:PUBLIC "Desktop"  # C:\Users\Public\Desktop
} else {
    $desktop = [Environment]::GetFolderPath("Desktop")
}
$lnkPath = Join-Path $desktop "Filipin Quant.lnk"

if (-not (Test-Path $target)) {
    Write-Host "ERROR: launch_filipin_gui.bat not found next to this script." -ForegroundColor Red
    Read-Host "Press Enter to exit"; exit 1
}

try {
    $ws = New-Object -ComObject WScript.Shell
    $sc = $ws.CreateShortcut($lnkPath)
    $sc.TargetPath       = $target
    $sc.WorkingDirectory = $workdir
    $sc.Description       = "Launch the Filipin quantification GUI"
    if (Test-Path $icon) { $sc.IconLocation = $icon }
    $sc.Save()
} catch {
    Write-Host "ERROR: could not write the shortcut to:" -ForegroundColor Red
    Write-Host "  $lnkPath"
    if ($AllUsers) { Write-Host "The Public Desktop needs admin rights - run this from an admin account." -ForegroundColor Yellow }
    Read-Host "Press Enter to exit"; exit 1
}

Write-Host ""
Write-Host "Created shortcut:" -ForegroundColor Green
Write-Host "  $lnkPath"
if ($AllUsers) { Write-Host "It will appear on the Desktop for ALL users (including standard users)." }
Write-Host "Double-click 'Filipin Quant' to start the GUI."
Write-Host ""
Read-Host "Press Enter to close"
