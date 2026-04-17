param(
    [string]$StudyName = "",
    [string]$PromptFile = "Prompt.md",
    [int]$PollSeconds = 20,
    [switch]$NoWatcher,
    [switch]$NoLiveView
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$Root = Split-Path -Parent $PSScriptRoot
$PromptPath = Join-Path $Root $PromptFile
$QuickstartScript = Join-Path (Join-Path $Root "tools") "research_loop.ps1"
$WatcherScript = Join-Path (Join-Path $Root "tools") "auto_loop.ps1"
$LiveViewScript = Join-Path (Join-Path $Root "tools") "agent_live_view.ps1"
$StudyDir = $null
$RunDir = $null
$StateFile = $null
$ReviewRequestFile = $null
$ReviewDirFile = $null

function Get-Slug {
    param([string]$Text)

    $value = $Text.ToLowerInvariant()
    $value = [regex]::Replace($value, '[^a-z0-9]+', '_')
    $value = [regex]::Replace($value, '_+', '_').Trim('_')
    if ($value.Length -gt 80) {
        $value = $value.Substring(0, 80).TrimEnd('_')
    }
    if ([string]::IsNullOrWhiteSpace($value)) {
        return 'auto_research_study'
    }
    return $value
}

function Resolve-PromptPath {
    param([string]$RequestedPromptFile)

    $requestedPath = Join-Path $Root $RequestedPromptFile
    if (Test-Path -LiteralPath $requestedPath) {
        return $requestedPath
    }

    if ($RequestedPromptFile -ieq 'Prompt.md') {
        $legacyPath = Join-Path $Root 'Auto_Research_Prompt.md'
        if (Test-Path -LiteralPath $legacyPath) {
            return $legacyPath
        }
    }

    return $requestedPath
}

function Resolve-StudyName {
    param(
        [string]$PromptText,
        [string]$ExplicitStudyName
    )

    if (-not [string]::IsNullOrWhiteSpace($ExplicitStudyName)) {
        return Get-Slug $ExplicitStudyName
    }

    $lines = @($PromptText -split "`r?`n")

    foreach ($line in $lines) {
        if ($line -match '^(StudyName|Study Name)\s*:\s*(.+?)\s*$') {
            return Get-Slug $Matches[2]
        }
    }

    foreach ($line in $lines) {
        $trimmed = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trimmed)) {
            continue
        }

        if ($trimmed -match '^#+\s*(.+?)\s*$') {
            $heading = $Matches[1].Trim()
            if ($heading -match '^[A-Za-z ]*Auto Research Prompt[A-Za-z ]*$') {
                continue
            }
            if ($heading -match '^[A-Za-z ]*AI Research Prompt\s*:\s*(.+)$') {
                return Get-Slug $Matches[1]
            }
            return Get-Slug $heading
        }
    }

    foreach ($line in $lines) {
        $trimmed = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trimmed)) {
            continue
        }
        return Get-Slug $trimmed
    }

    return 'auto_research_study'
}

function Quote-PSValue {
    param([string]$Text)

    return "'" + $Text.Replace("'", "''") + "'"
}

function Read-StudyStatus {
    if (-not $StateFile -or -not (Test-Path -LiteralPath $StateFile)) {
        return $null
    }

    try {
        $state = Get-Content -LiteralPath $StateFile -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($state -and $state.status) {
            return ([string]$state.status).ToUpperInvariant()
        }
    } catch {
    }

    return $null
}

function Read-ReviewDecision {
    if (-not $ReviewDirFile -or -not (Test-Path -LiteralPath $ReviewDirFile)) {
        return $null
    }

    try {
        $content = Get-Content -LiteralPath $ReviewDirFile -Raw -Encoding UTF8
    } catch {
        return $null
    }

    if ($content -match '##\s*Decision\s*[\r\n]+\s*(APPROVE|REVISE|NEEDS_REWORK)') {
        return $Matches[1]
    }
    if ($content -match '\bAPPROVE\b') { return 'APPROVE' }
    if ($content -match '\bNEEDS_REWORK\b') { return 'NEEDS_REWORK' }
    if ($content -match '\bREVISE\b') { return 'REVISE' }
    return $null
}

function Get-InitialStagePrompt {
    param([string]$ResolvedStudyName)

    $relStudy = "studies/$ResolvedStudyName"
    $relRun = "task_runs/$ResolvedStudyName"
    $status = Read-StudyStatus
    $decision = Read-ReviewDecision

    if ($decision -eq 'APPROVE' -or $status -eq 'APPROVED') {
        return "@execution-engineer study=$relStudy run=$relRun phase=polish"
    }

    if ($decision -eq 'REVISE' -or $decision -eq 'NEEDS_REWORK' -or $status -eq 'REVISION_REQUESTED') {
        return "@execution-engineer study=$relStudy run=$relRun phase=implement"
    }

    if ((Test-Path -LiteralPath $ReviewRequestFile) -or $status -eq 'REVIEW_REQUESTED' -or $status -eq 'READY_FOR_REVIEW') {
        return "@science-director study=$relStudy run=$relRun phase=review"
    }

    switch ($status) {
        'PLANNED' {
            return "@execution-engineer study=$relStudy run=$relRun phase=implement"
        }
        'IMPLEMENTED' {
            return "@execution-engineer study=$relStudy run=$relRun phase=validate"
        }
        'REPORTING' {
            return "@execution-engineer study=$relStudy run=$relRun phase=report"
        }
        'COMPLETE' {
            return $null
        }
        default {
            return "@execution-engineer study=$relStudy run=$relRun phase=plan"
        }
    }
}

function Write-QuickstartSummary {
    param([object[]]$Lines)

    $skipPromptBlock = $false
    foreach ($lineObject in @($Lines)) {
        $line = [string]$lineObject

        if ($skipPromptBlock) {
            if ($line -match '^Model to select:') {
                $skipPromptBlock = $false
            } else {
                continue
            }
        }

        if ($line -match 'PASTE THIS INTO COPILOT CHAT') {
            $skipPromptBlock = $true
            continue
        }

        if ($line -match '^\s*Clipboard: agent invocation copied\.$') {
            continue
        }

        if ($line -match '^Command is in your clipboard\. Ctrl\+V into Copilot Chat\.$') {
            continue
        }

        if ($line -match '^Model to select:') {
            continue
        }

        if ($line -match '^Context file:') {
            continue
        }

        Write-Host $line
    }
}

function Get-ScriptProcesses {
    param(
        [string]$ScriptPath,
        [string]$StudyNameFilter
    )

    $scriptPattern = [regex]::Escape((Split-Path -Leaf $ScriptPath))
    $studyPattern = [regex]::Escape($StudyNameFilter)

    return @(
        Get-CimInstance Win32_Process |
            Where-Object {
                $_.Name -eq 'powershell.exe' -and
                $_.CommandLine -and
                $_.CommandLine -match $scriptPattern -and
                $_.CommandLine -match $studyPattern
            }
    )
}

function Stop-ScriptProcesses {
    param(
        [string]$Label,
        [string]$ScriptPath,
        [string]$StudyNameFilter
    )

    $existing = Get-ScriptProcesses -ScriptPath $ScriptPath -StudyNameFilter $StudyNameFilter
    foreach ($proc in @($existing)) {
        try {
            Stop-Process -Id $proc.ProcessId -Force -ErrorAction Stop
            Write-Host (("Stopped existing " + $Label + " PID ") + $proc.ProcessId) -ForegroundColor DarkYellow
        } catch {
            Write-Host (("Could not stop existing " + $Label + " PID ") + $proc.ProcessId + ": " + $_.Exception.Message) -ForegroundColor Yellow
        }
    }
}

function Start-VisibleScriptWindow {
    param(
        [string]$Label,
        [string]$WindowTitle,
        [string]$ScriptPath,
        [string[]]$ScriptArguments
    )

    $scriptCall = "& powershell.exe -NoProfile -ExecutionPolicy Bypass -File " + (Quote-PSValue $ScriptPath)
    foreach ($token in @($ScriptArguments)) {
        if ($token.StartsWith("-")) {
            $scriptCall += " " + $token
        } else {
            $scriptCall += " " + (Quote-PSValue $token)
        }
    }

    $command = @"
`$host.UI.RawUI.WindowTitle = $(Quote-PSValue $WindowTitle)
try {
    $scriptCall
} catch {
    Write-Host ('Unhandled launcher error: ' + `$_.Exception.Message) -ForegroundColor Red
}
`$exitCodeVar = Get-Variable -Name LASTEXITCODE -ErrorAction SilentlyContinue
`$exitCode = if (`$exitCodeVar) { [int]`$exitCodeVar.Value } else { 0 }
if (`$exitCode -ne 0) {
    Write-Host ''
    Write-Host ($(Quote-PSValue ($Label + ' exited with code ')) + `$exitCode) -ForegroundColor Red
    Read-Host 'Press Enter to close' | Out-Null
}
"@

    return Start-Process -FilePath "powershell.exe" -ArgumentList @(
        "-NoExit",
        "-Command",
        $command
    ) -WorkingDirectory $Root -WindowStyle Normal -PassThru
}

$PromptPath = Resolve-PromptPath -RequestedPromptFile $PromptFile

Write-Host "" 
Write-Host "=== Auto Agent Study Launcher ===" -ForegroundColor Cyan
Write-Host "Prompt    : $PromptFile" -ForegroundColor White
Write-Host "" 

if (-not (Test-Path -LiteralPath $PromptPath)) {
    Write-Host "Prompt file not found: $PromptPath" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $QuickstartScript)) {
    Write-Host "Missing script: $QuickstartScript" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $WatcherScript)) {
    Write-Host "Missing script: $WatcherScript" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path -LiteralPath $LiveViewScript)) {
    Write-Host "Missing script: $LiveViewScript" -ForegroundColor Red
    exit 1
}

$studyGoal = Get-Content -LiteralPath $PromptPath -Raw -Encoding UTF8
if ([string]::IsNullOrWhiteSpace($studyGoal)) {
    Write-Host "Prompt file is empty: $PromptPath" -ForegroundColor Red
    exit 1
}

$StudyName = Resolve-StudyName -PromptText $studyGoal -ExplicitStudyName $StudyName
$StudyDir = Join-Path (Join-Path $Root "studies") $StudyName
$RunDir = Join-Path (Join-Path $Root "task_runs") $StudyName
$StateFile = Join-Path $StudyDir "study_state.json"
$ReviewRequestFile = Join-Path $RunDir "REVIEW_REQUEST.md"
$ReviewDirFile = Join-Path $RunDir "REVIEW_DIRECTIVE.md"
Write-Host "StudyName : $StudyName" -ForegroundColor White
Write-Host "PromptPath : $PromptPath" -ForegroundColor White
Write-Host "" 

Write-Host "Running quickstart (init or extend study)..." -ForegroundColor Yellow
$LASTEXITCODE = 0
$quickstartOutput = @(& $QuickstartScript -Action quickstart -StudyName $StudyName -StudyGoal $studyGoal 2>&1)
$quickstartExitCode = if ($LASTEXITCODE) { [int]$LASTEXITCODE } else { 0 }
Write-QuickstartSummary -Lines $quickstartOutput
if ($quickstartExitCode -ne 0 -or -not $?) {
    Write-Host "Quickstart failed. Fix errors above and retry." -ForegroundColor Red
    exit 1
}

$initialStagePrompt = Get-InitialStagePrompt -ResolvedStudyName $StudyName
if (-not [string]::IsNullOrWhiteSpace($initialStagePrompt)) {
    $copiedInitialPrompt = $false
    try {
        $initialStagePrompt | Set-Clipboard
        $copiedInitialPrompt = $true
    } catch {
    }

    Write-Host "" 
    Write-Host "Staged workflow prompt:" -ForegroundColor Yellow
    Write-Host ("  " + $initialStagePrompt) -ForegroundColor White
    if ($copiedInitialPrompt) {
        Write-Host "Copied the staged prompt to the clipboard." -ForegroundColor Green
    }
    Write-Host "This launcher uses the staged watcher flow; ignore any one-shot @research-loop prompt from quickstart." -ForegroundColor DarkYellow
}

if ($NoWatcher) {
    Write-Host "Watcher launch skipped due to -NoWatcher." -ForegroundColor Yellow
} else {
    Stop-ScriptProcesses -Label "watcher" -ScriptPath $WatcherScript -StudyNameFilter $StudyName
    Write-Host "Starting watcher in a new PowerShell window..." -ForegroundColor Yellow
    $watcherProcess = Start-VisibleScriptWindow -Label "Watcher" -WindowTitle ("cQED Watcher - " + $StudyName) -ScriptPath $WatcherScript -ScriptArguments @(
        "-StudyName",
        $StudyName,
        "-PollSeconds",
        $PollSeconds.ToString()
    )
    Write-Host ("Watcher PID: " + $watcherProcess.Id) -ForegroundColor DarkGray
}

if ($NoLiveView) {
    Write-Host "Live view launch skipped due to -NoLiveView." -ForegroundColor Yellow
} else {
    Stop-ScriptProcesses -Label "live-view" -ScriptPath $LiveViewScript -StudyNameFilter $StudyName
    Write-Host "Starting live-view monitor in a new PowerShell window..." -ForegroundColor Yellow
    $liveViewProcess = Start-VisibleScriptWindow -Label "Live-view" -WindowTitle ("cQED Live View - " + $StudyName) -ScriptPath $LiveViewScript -ScriptArguments @(
        "-StudyName",
        $StudyName
    )
    Write-Host ("Live-view PID: " + $liveViewProcess.Id) -ForegroundColor DarkGray
}

Write-Host "" 
Write-Host "Next steps:" -ForegroundColor Green
if ($NoWatcher) {
    Write-Host "1) Paste the staged workflow prompt above into Copilot Chat with the requested model." -ForegroundColor Green
} else {
    Write-Host "1) Keep the watcher window open; it will mirror the same staged prompt and advance the loop." -ForegroundColor Green
}
if ($NoLiveView) {
    Write-Host "2) Keep GitHub Copilot Chat open to watch the agent stream in real time." -ForegroundColor Green
} else {
    Write-Host "2) Keep the live-view window open for file-level agent output updates." -ForegroundColor Green
}
Write-Host "3) Keep GitHub Copilot Chat open to see the agent stream in real time." -ForegroundColor Green
Write-Host "4) Use the staged prompt sequence (plan -> implement -> validate -> report -> review -> polish)." -ForegroundColor Green
Write-Host "" 
Write-Host "Done." -ForegroundColor Green
