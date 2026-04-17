@echo off
setlocal

REM One-click launcher for long research prompts.
REM Starts quickstart, watcher, and live-view monitor windows.
REM Study name is auto-derived from Prompt.md or Auto_Research_Prompt.md.

set "ROOT=%~dp0"
set "PROMPT_FILE=Prompt.md"
set "POLL_SECONDS=20"

if not exist "%ROOT%%PROMPT_FILE%" (
  set "PROMPT_FILE=Auto_Research_Prompt.md"
)

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%ROOT%tools\auto_agent_study.ps1" -PromptFile "%PROMPT_FILE%" -PollSeconds %POLL_SECONDS%
if errorlevel 1 (
  echo.
  echo Launcher failed. See errors above.
  pause
  exit /b 1
)

echo.
echo Launcher completed.
pause
exit /b 0
