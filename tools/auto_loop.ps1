<#
.SYNOPSIS
    Background watcher that detects research loop state changes and prompts the user to paste into Copilot Chat.

.DESCRIPTION
    Polls task_runs/<study>/ for state-file signals. When a stage completes, the script:
    1. Generates the appropriate agent prompt
    2. Copies it to clipboard
    3. Shows instructions for which model to select in Copilot Chat
    4. Waits for the next stage output file

    All agents are run through GitHub Copilot Chat. Users:
    - Open GitHub Copilot Chat in VS Code
    - Select the appropriate model (Opus 4.6 or Codex 5.4 xHigh)
    - Paste the prompt (auto-copied to clipboard)
    - Wait for results
    - Optional: save output to the state file if needed

    Stages and their models:
      Stage 1 -- Execute:    Opus 4.6 (via Copilot Chat)
      Stage 2 -- Review:     Codex 5.4 xHigh (via Copilot Chat)
      Stage 3 -- Refine:     Opus 4.6 (via Copilot Chat)
      Stage 4 -- Polish:     Opus 4.6 (via Copilot Chat)

.PARAMETER StudyName
    Study folder name under studies/  (e.g. chi_sweep)

.PARAMETER PollSeconds
    How often (in seconds) to check for new signals. Default: 20.
    Increase if your Box sync is slow.

.PARAMETER Model
    Claude model to use for execution/polish. Default: claude-opus-4-6.

.PARAMETER DryRun
    Print what would be run but do not actually call the claude CLI.

.EXAMPLE
    # Start the watcher for a study
    powershell -ExecutionPolicy Bypass -File tools\auto_loop.ps1 -StudyName chi_sweep

    # Dry-run to see the prompts without firing claude
    powershell -ExecutionPolicy Bypass -File tools\auto_loop.ps1 -StudyName chi_sweep -DryRun
#>
param(
    [Parameter(Mandatory=$true)]
    [string]$StudyName,

    [int]$PollSeconds = 20,

    [string]$Model = "claude-opus-4-6",

    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
$Root           = Split-Path -Parent $PSScriptRoot
$StudyDir       = Join-Path (Join-Path $Root "studies") $StudyName
$RunDir         = Join-Path (Join-Path $Root "task_runs") $StudyName
$ConfigFile     = Join-Path $Root "research_config.json"
$AgentFile      = Join-Path (Join-Path (Join-Path $Root ".github") "agents") "execution-engineer.agent.md"

$ScienceDirective = Join-Path $RunDir "SCIENCE_DIRECTIVE.md"
$ExecutionSummary = Join-Path $RunDir "EXECUTION_SUMMARY.md"
$ReviewRequestFile = Join-Path $RunDir "REVIEW_REQUEST.md"
$ReviewDirFile     = Join-Path $RunDir "REVIEW_DIRECTIVE.md"
$FollowupFile      = Join-Path $RunDir "FOLLOWUP_PROMPT.md"
$PolishDoneFile    = Join-Path $RunDir "POLISH_COMPLETE.md"
$StateFile         = Join-Path $StudyDir "study_state.json"

$RelStudy       = "studies/$StudyName"
$RelRun         = "task_runs/$StudyName"

# ---------------------------------------------------------------------------
# Load executor model from research_config.json (overrides -Model if present)
# ---------------------------------------------------------------------------
if (Test-Path $ConfigFile) {
    try {
        $cfg = Get-Content $ConfigFile -Raw | ConvertFrom-Json
        if ($cfg.models.execution.model_id -and $cfg.models.execution.model_id -notlike "*Codex*") {
            $Model = $cfg.models.execution.model_id
        }
    } catch { }
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
function Write-Banner {
    param([string]$Text, [string]$Color = "Cyan")
    $line = "=" * 70
    Write-Host ""
    Write-Host $line -ForegroundColor $Color
    Write-Host "  $Text" -ForegroundColor $Color
    Write-Host $line -ForegroundColor $Color
    Write-Host ""
}

function Write-Step {
    param([string]$Text, [string]$Color = "White")
    $ts = (Get-Date).ToString("HH:mm:ss")
    Write-Host "[$ts] $Text" -ForegroundColor $Color
}

function Invoke-UserAlert {
    param([string]$Message)
    # Beep
    [System.Console]::Beep(800, 300)
    Start-Sleep -Milliseconds 100
    [System.Console]::Beep(1000, 300)

    # Windows toast notification (best-effort)
    try {
        $null = [Windows.UI.Notifications.ToastNotificationManager, Windows.UI.Notifications, ContentType=WindowsRuntime]
        $template = [Windows.UI.Notifications.ToastNotificationManager]::GetTemplateContent(
            [Windows.UI.Notifications.ToastTemplateType]::ToastText02)
        $template.GetElementsByTagName("text")[0].InnerText = "cQED Research Loop"
        $template.GetElementsByTagName("text")[1].InnerText = $Message
        $notifier = [Windows.UI.Notifications.ToastNotificationManager]::CreateToastNotifier("ResearchLoop")
        $notifier.Show([Windows.UI.Notifications.ToastNotification]::new($template))
    } catch { } # Silently skip if toast not available
}

function Test-CLIs {
    # No CLI verification needed — all work goes through Copilot Chat UI
}

function Copy-ToClipboard {
    param([string]$Text)
    try {
        $Text | Set-Clipboard
        return $true
    } catch {
        return $false
    }
}

function Get-AgentInstructions {
    # Read the science-director agent file, strip YAML frontmatter
    if (-not (Test-Path $AgentFile)) {
        return "You are the Critical Reviewer for a cQED research loop. Read AGENTS.md for full instructions."
    }
    $raw = Get-Content $AgentFile -Raw -Encoding UTF8
    # Strip YAML frontmatter (--- ... ---)
    $raw = [regex]::Replace($raw, '(?s)^---.*?---\s*', '')
    return $raw.Trim()
}

function Get-FileTimestamp {
    param([string]$Path)
    if (Test-Path $Path) {
        return (Get-Item $Path).LastWriteTimeUtc
    }
    return [datetime]::MinValue
}

function Read-Decision {
    # Parse REVIEW_DIRECTIVE.md for the decision line
    if (-not (Test-Path $ReviewDirFile)) { return $null }
    $content = Get-Content $ReviewDirFile -Raw -Encoding UTF8
    if ($content -match '##\s*Decision\s*[\r\n]+\s*(APPROVE|REVISE|NEEDS_REWORK)') {
        return $Matches[1]
    }
    # Also check bare decision keyword anywhere in the file (fallback)
    if ($content -match '\bAPPROVE\b') { return 'APPROVE' }
    if ($content -match '\bNEEDS_REWORK\b') { return 'NEEDS_REWORK' }
    if ($content -match '\bREVISE\b') { return 'REVISE' }
    return $null
}

function Read-StudyStatus {
    if (-not (Test-Path $StateFile)) { return $null }
    try {
        $state = Get-Content $StateFile -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($state -and $state.status) {
            return ([string]$state.status).ToUpperInvariant()
        }
    } catch {
    }
    return $null
}

function New-LoopAction {
    param(
        [string]$Key,
        [string]$Phase,
        [string]$WaitingFor,
        [string]$Message,
        [string]$Color = 'Yellow'
    )

    return [PSCustomObject]@{
        Key = $Key
        Phase = $Phase
        WaitingFor = $WaitingFor
        Message = $Message
        Color = $Color
    }
}

function Get-NextLoopAction {
    $status = Read-StudyStatus
    $decision = Read-Decision
    $hasScienceDirective = Test-Path $ScienceDirective
    $hasExecutionSummary = Test-Path $ExecutionSummary
    $hasReviewRequest = Test-Path $ReviewRequestFile
    $hasFollowup = Test-Path $FollowupFile

    if ((Test-Path $PolishDoneFile) -or $status -eq 'COMPLETE') {
        return New-LoopAction -Key 'complete' -Phase '' -WaitingFor 'none' -Message 'Study is complete.' -Color 'Green'
    }

    switch ($status) {
        'APPROVED' {
            return New-LoopAction -Key 'polish' -Phase 'polish' -WaitingFor 'POLISH_COMPLETE.md' -Message 'Review approved. Showing polish prompt...' -Color 'Green'
        }
        'REVIEW_REQUESTED' {
            return New-LoopAction -Key 'review' -Phase 'review' -WaitingFor 'REVIEW_DIRECTIVE.md' -Message 'REVIEW_REQUEST.md found. Showing Codex review prompt...' -Color 'Yellow'
        }
        'READY_FOR_REVIEW' {
            return New-LoopAction -Key 'review' -Phase 'review' -WaitingFor 'REVIEW_DIRECTIVE.md' -Message 'Study is ready for review. Showing Codex review prompt...' -Color 'Yellow'
        }
        'REVIEWING' {
            return New-LoopAction -Key 'review' -Phase 'review' -WaitingFor 'REVIEW_DIRECTIVE.md' -Message 'Review is in progress. Waiting for REVIEW_DIRECTIVE.md...' -Color 'Yellow'
        }
        'REVISION_REQUESTED' {
            return New-LoopAction -Key 'refine' -Phase 'implement' -WaitingFor 'updated EXECUTION_SUMMARY.md and REVIEW_REQUEST.md' -Message 'Reviewer requested another implementation pass.' -Color 'Yellow'
        }
        'INITIALIZED' {
            if ($hasScienceDirective) {
                return New-LoopAction -Key 'implement' -Phase 'implement' -WaitingFor 'EXECUTION_SUMMARY.md' -Message 'SCIENCE_DIRECTIVE.md exists. Showing implementation prompt...' -Color 'Yellow'
            }
            return New-LoopAction -Key 'plan' -Phase 'plan' -WaitingFor 'SCIENCE_DIRECTIVE.md' -Message 'Study initialized. Showing planning prompt...' -Color 'Yellow'
        }
        'PLANNING' {
            if ($hasScienceDirective) {
                return New-LoopAction -Key 'implement' -Phase 'implement' -WaitingFor 'EXECUTION_SUMMARY.md' -Message 'SCIENCE_DIRECTIVE.md detected. Showing implementation prompt...' -Color 'Yellow'
            }
            return New-LoopAction -Key 'plan' -Phase 'plan' -WaitingFor 'SCIENCE_DIRECTIVE.md' -Message 'Planning is in progress. Waiting for SCIENCE_DIRECTIVE.md...' -Color 'Yellow'
        }
        'PLANNED' {
            return New-LoopAction -Key 'implement' -Phase 'implement' -WaitingFor 'EXECUTION_SUMMARY.md' -Message 'Study is planned. Showing implementation prompt...' -Color 'Yellow'
        }
        'IMPLEMENTING' {
            return New-LoopAction -Key 'implement' -Phase 'implement' -WaitingFor 'EXECUTION_SUMMARY.md' -Message 'Implementation is in progress. Waiting for EXECUTION_SUMMARY.md...' -Color 'Yellow'
        }
        'IMPLEMENTED' {
            return New-LoopAction -Key 'validate' -Phase 'validate' -WaitingFor 'study_state.json -> REPORTING' -Message 'Implementation complete. Showing validation prompt...' -Color 'Yellow'
        }
        'VALIDATING' {
            return New-LoopAction -Key 'validate' -Phase 'validate' -WaitingFor 'study_state.json -> REPORTING' -Message 'Validation is in progress. Waiting for report-phase readiness...' -Color 'Yellow'
        }
        'REPORTING' {
            return New-LoopAction -Key 'report' -Phase 'report' -WaitingFor 'REVIEW_REQUEST.md' -Message 'Validation is complete. Showing report prompt...' -Color 'Yellow'
        }
    }

    if ($decision -eq 'APPROVE') {
        return New-LoopAction -Key 'polish' -Phase 'polish' -WaitingFor 'POLISH_COMPLETE.md' -Message 'Review approved. Showing polish prompt...' -Color 'Green'
    }

    if ($decision -eq 'REVISE' -or $decision -eq 'NEEDS_REWORK') {
        return New-LoopAction -Key 'refine' -Phase 'implement' -WaitingFor 'updated EXECUTION_SUMMARY.md and REVIEW_REQUEST.md' -Message 'Reviewer requested another implementation pass.' -Color 'Yellow'
    }

    if ($hasReviewRequest) {
        return New-LoopAction -Key 'review' -Phase 'review' -WaitingFor 'REVIEW_DIRECTIVE.md' -Message 'Review handoff exists. Showing review prompt...' -Color 'Yellow'
    }

    if ($hasFollowup -and -not $hasReviewRequest) {
        return New-LoopAction -Key 'refine' -Phase 'implement' -WaitingFor 'updated EXECUTION_SUMMARY.md and REVIEW_REQUEST.md' -Message 'FOLLOWUP_PROMPT.md found. Showing refinement prompt...' -Color 'Yellow'
    }

    if (-not $hasScienceDirective) {
        return New-LoopAction -Key 'plan' -Phase 'plan' -WaitingFor 'SCIENCE_DIRECTIVE.md' -Message 'No SCIENCE_DIRECTIVE.md found. Showing planning prompt...' -Color 'Yellow'
    }

    if (-not $hasExecutionSummary) {
        return New-LoopAction -Key 'implement' -Phase 'implement' -WaitingFor 'EXECUTION_SUMMARY.md' -Message 'No EXECUTION_SUMMARY.md found. Showing implementation prompt...' -Color 'Yellow'
    }

    if (-not $hasReviewRequest) {
        return New-LoopAction -Key 'report' -Phase 'report' -WaitingFor 'REVIEW_REQUEST.md' -Message 'Execution summary exists but review handoff is missing. Showing report prompt...' -Color 'Yellow'
    }

    return New-LoopAction -Key 'review' -Phase 'review' -WaitingFor 'REVIEW_DIRECTIVE.md' -Message 'Review handoff exists. Showing review prompt...' -Color 'Yellow'
}

# ---------------------------------------------------------------------------
# Stage prompt display
# ---------------------------------------------------------------------------

function Show-CopilotPrompt {
    param(
        [string]$Prompt,
        [string]$Stage,
        [string]$ModelToSelect
    )
    
    $copied = $false
    if (-not $DryRun) {
        $copied = Copy-ToClipboard $Prompt
    }
    Invoke-UserAlert "Action needed: paste Copilot prompt for $Stage"

    Write-Banner "ACTION REQUIRED -- PASTE INTO GITHUB COPILOT CHAT" "Magenta"
    Write-Host "  Stage: $Stage" -ForegroundColor Magenta
    Write-Host "  Model to select: $ModelToSelect" -ForegroundColor Magenta
    Write-Host ""
    
    if ($copied) {
        Write-Host "  [YES] Prompt copied to clipboard. Ctrl+V in Copilot Chat." -ForegroundColor Green
    } elseif ($DryRun) {
        Write-Host "  [DRY RUN] Prompt not copied because -DryRun was requested." -ForegroundColor Yellow
    } else {
        Write-Host "  [!] Clipboard copy failed. Prompt shown below:" -ForegroundColor Yellow
    }
    Write-Host ""
    Write-Host "  +-------------------------------------------------------------+" -ForegroundColor White
    foreach ($line in $Prompt -split "`n") {
        Write-Host ("  | " + $line) -ForegroundColor White
    }
    Write-Host "  +-------------------------------------------------------------+" -ForegroundColor White
    Write-Host ""
    Write-Host "  Waiting for output file to appear..." -ForegroundColor Gray
}

# ---------------------------------------------------------------------------
# Stage runners (generate prompts and show them)
# ---------------------------------------------------------------------------

function Invoke-Execute {
    param([string]$Phase = "plan")
    
    $relStudy2 = "studies/$StudyName"
    $relRun2   = "task_runs/$StudyName"
    $prompt = "@execution-engineer study=$relStudy2 run=$relRun2 phase=$Phase"
    
    Show-CopilotPrompt -Prompt $prompt -Stage "Stage 1 -- Execute ($Phase)" -ModelToSelect "Opus 4.6"
}

function Invoke-Review {
    $relStudy2 = "studies/$StudyName"
    $relRun2   = "task_runs/$StudyName"
    $prompt = "@science-director study=$relStudy2 run=$relRun2 phase=review"
    
    Show-CopilotPrompt -Prompt $prompt -Stage "Stage 2 -- Review" -ModelToSelect "Codex 5.4 xHigh"
}

function Invoke-Refine {
    $relStudy2 = "studies/$StudyName"
    $relRun2   = "task_runs/$StudyName"
    $prompt = "@execution-engineer study=$relStudy2 run=$relRun2 phase=implement"
    
    if (Test-Path $FollowupFile) {
        $followup = Get-Content $FollowupFile -Raw -Encoding UTF8
        $prompt += "`n`n$followup"
    }
    
    Show-CopilotPrompt -Prompt $prompt -Stage "Stage 3 -- Refine" -ModelToSelect "Opus 4.6"
}

function Invoke-Polish {
    $relStudy2 = "studies/$StudyName"
    $relRun2   = "task_runs/$StudyName"
    $prompt = "@execution-engineer study=$relStudy2 run=$relRun2 phase=polish"
    
    Show-CopilotPrompt -Prompt $prompt -Stage "Stage 4 -- Polish" -ModelToSelect "Opus 4.6"
}

function Show-LoopAction {
    param([PSObject]$Action)

    switch ($Action.Key) {
        'plan' {
            Invoke-Execute -Phase 'plan'
        }
        'implement' {
            Invoke-Execute -Phase 'implement'
        }
        'validate' {
            Invoke-Execute -Phase 'validate'
        }
        'report' {
            Invoke-Execute -Phase 'report'
        }
        'review' {
            Invoke-Review
        }
        'refine' {
            Invoke-Refine
        }
        'polish' {
            Invoke-Polish
        }
    }
}

# ---------------------------------------------------------------------------
# Main watch loop
# ---------------------------------------------------------------------------

Write-Banner "cQED RESEARCH LOOP -- AUTO WATCHER" "Cyan"
Write-Step "Study:        $StudyName" "White"
Write-Step "Study dir:    $RelStudy" "White"
Write-Step "Run dir:      $RelRun" "White"
Write-Step "Exec model:   $Model (Opus 4.6 via Copilot Chat)" "White"
Write-Step "Review model: Codex 5.4 xHigh (via Copilot Chat)" "White"
Write-Step "Poll interval: ${PollSeconds}s" "White"
Write-Step "All prompts will be copied to clipboard and require manual paste into Copilot Chat." "Gray"
if ($DryRun) { Write-Step "MODE: DRY RUN (prompts will be shown but not copied)" "DarkYellow" }
Write-Host ""

# Verify CLI availability
Test-CLIs

# Validate paths
if (-not (Test-Path $StudyDir)) {
    Write-Host "Study directory not found: $StudyDir" -ForegroundColor Red
    Write-Host "Run 'Research: New Study' VS Code task first, or use -Action quickstart." -ForegroundColor Red
    exit 1
}

# Ensure run dir exists
New-Item -ItemType Directory -Path $RunDir -Force | Out-Null

# --- Timestamps of last-processed files ---
$lastReviewRequestTime = Get-FileTimestamp $ReviewRequestFile
$lastReviewDirTime     = Get-FileTimestamp $ReviewDirFile
$lastStatus            = Read-StudyStatus
$currentAction         = Get-NextLoopAction
$lastActionKey         = $currentAction.Key

if ($currentAction.Key -eq 'complete') {
    Write-Banner "STUDY IS COMPLETE -- POLISH DONE" "Green"
    Write-Step "See: $RelRun/POLISH_COMPLETE.md" "Green"
    Write-Step "Report: $RelStudy/report/report.pdf" "Green"
    exit 0
}

Write-Step $currentAction.Message $currentAction.Color
Show-LoopAction -Action $currentAction

Write-Step "Watcher active. Ctrl+C to stop." "Gray"
Write-Host ""

# Track study directory activity for heartbeat
$StudySubdirs = @('scripts', 'data', 'figures', 'report', 'artifacts')
$lastStudyFileCount = 0
try {
    foreach ($sub in $StudySubdirs) {
        $subPath = Join-Path $StudyDir $sub
        if (Test-Path -LiteralPath $subPath) {
            $lastStudyFileCount += @(Get-ChildItem -LiteralPath $subPath -Recurse -File -ErrorAction SilentlyContinue).Count
        }
    }
} catch { }
$HeartbeatCycles = [math]::Max(1, [math]::Floor(300 / $PollSeconds))  # ~5min between heartbeats
$cycleCounter = 0

while ($true) {
    $now = Get-Date
    $cycleCounter++

    $currentReviewRequestTime = Get-FileTimestamp $ReviewRequestFile
    $currentReviewDirTime = Get-FileTimestamp $ReviewDirFile
    $currentStatus = Read-StudyStatus
    $currentAction = Get-NextLoopAction

    if ($currentStatus -ne $lastStatus) {
        $fromStatus = if ($lastStatus) { $lastStatus } else { '(none)' }
        $toStatus = if ($currentStatus) { $currentStatus } else { '(none)' }
        Write-Step ("Study status changed: " + $fromStatus + " -> " + $toStatus) "DarkCyan"
        $lastStatus = $currentStatus
    }

    if ($currentAction.Key -eq 'complete') {
        Write-Host ""
        Write-Banner "STUDY IS COMPLETE -- POLISH DONE" "Green"
        Write-Step "See: $RelRun/POLISH_COMPLETE.md" "Green"
        Write-Step "Report: $RelStudy/report/report.pdf" "Green"
        exit 0
    }

    $shouldPrompt = $false
    if ($currentAction.Key -ne $lastActionKey) {
        $shouldPrompt = $true
    } elseif ($currentAction.Key -eq 'review' -and $currentReviewRequestTime -gt $lastReviewRequestTime) {
        $shouldPrompt = $true
    } elseif (($currentAction.Key -eq 'refine' -or $currentAction.Key -eq 'polish') -and $currentReviewDirTime -gt $lastReviewDirTime) {
        $shouldPrompt = $true
    }

    if ($shouldPrompt) {
        Write-Step $currentAction.Message $currentAction.Color
        Show-LoopAction -Action $currentAction
        $lastActionKey = $currentAction.Key
    }

    $lastReviewRequestTime = $currentReviewRequestTime
    $lastReviewDirTime = $currentReviewDirTime
    # ---- Heartbeat with study activity info ----
    if (($cycleCounter % $HeartbeatCycles) -eq 0) {
        $currentFileCount = 0
        $recentFiles = @()
        try {
            foreach ($sub in $StudySubdirs) {
                $subPath = Join-Path $StudyDir $sub
                if (Test-Path -LiteralPath $subPath) {
                    $subFiles = @(Get-ChildItem -LiteralPath $subPath -Recurse -File -ErrorAction SilentlyContinue)
                    $currentFileCount += @($subFiles).Count
                    foreach ($f in @($subFiles)) {
                        if (($now - $f.LastWriteTime).TotalSeconds -lt ($HeartbeatCycles * $PollSeconds + 5)) {
                            $relName = $f.FullName
                            if ($f.FullName.StartsWith($StudyDir)) {
                                $relName = $f.FullName.Substring($StudyDir.Length).TrimStart('\','/')
                            }
                            $recentFiles += $relName
                        }
                    }
                }
            }
        } catch { }

        $ts = $now.ToString("HH:mm:ss")
        $waitingFor = $currentAction.WaitingFor

        $newFilesDelta = $currentFileCount - $lastStudyFileCount
        $deltaStr = ""
        if ($newFilesDelta -gt 0) {
            $deltaStr = ", +" + $newFilesDelta + " new files"
        }

        Write-Host ("[$ts] Watching... ($currentFileCount study files" + $deltaStr + ") waiting for: $waitingFor") -ForegroundColor DarkGray

        if (@($recentFiles).Count -gt 0) {
            foreach ($rf in @($recentFiles) | Select-Object -First 5) {
                Write-Host ("       ~ $rf") -ForegroundColor DarkCyan
            }
            if (@($recentFiles).Count -gt 5) {
                Write-Host ("       ... and " + (@($recentFiles).Count - 5) + " more") -ForegroundColor DarkCyan
            }
        }

        $lastStudyFileCount = $currentFileCount
    }

    Start-Sleep -Seconds $PollSeconds
}
