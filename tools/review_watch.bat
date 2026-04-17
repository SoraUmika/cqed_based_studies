@echo off
REM ============================================================
REM  cQED Research Loop — Auto-Watcher Launcher
REM  Usage: double-click this file, or run from terminal:
REM    tools\review_watch.bat <study_name> [poll_seconds]
REM
REM  This starts the background watcher that automatically
REM  runs the Opus 4.6 execution, refinement, and polish stages
REM  via the claude CLI. Only the Codex review stage (Stage 2)
REM  requires manual action (paste prompt into Copilot).
REM
REM  When human action is needed (paste into Copilot), the
REM  watcher beeps, shows the prompt, and copies it to clipboard.
REM ============================================================

setlocal

SET STUDY=%1
SET POLL=%2

IF "%STUDY%"=="" (
    echo.
    echo  Usage: review_watch.bat ^<study_name^> [poll_seconds]
    echo.
    echo  Example: review_watch.bat chi_sweep
    echo           review_watch.bat chi_sweep 30
    echo.
    SET /P STUDY=Enter study name (folder under studies/):
)

IF "%POLL%"=="" SET POLL=20

echo.
echo  Starting auto-watcher for study: %STUDY%
echo  Poll interval: %POLL% seconds
echo  Press Ctrl+C to stop.
echo.

REM Change to the workspace root (parent of this script's directory)
cd /d "%~dp0.."

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0auto_loop.ps1" -StudyName "%STUDY%" -PollSeconds %POLL%

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo  Watcher exited with error %ERRORLEVEL%.
    pause
)
