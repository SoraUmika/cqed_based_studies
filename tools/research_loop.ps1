<#
.SYNOPSIS
    Orchestrates the continuous two-model research loop for cQED studies.

.DESCRIPTION
    Manages the lifecycle of a research study through a Science Director (planning/review)
    and Execution Engineer (implementation/reporting) loop. Handles bootstrapping, state
    management, phase transitions, failure recovery, and iteration limits.

    Reads research_config.json at the workspace root for all configurable settings.

.PARAMETER Action
    quickstart     — ONE COMMAND: init (if needed) + build the full-loop prompt + copy to clipboard
    init           — Bootstrap a new study
    status         — Show current phase and open tasks
    plan           — Transition to PLANNING, print Science Director invocation
    execute        — Transition to IMPLEMENTING, print Execution Engineer invocation
    review         — Transition to REVIEWING, print Science Director review invocation
    validate       — Transition to VALIDATING, print Execution Engineer invocation
    report         — Transition to REPORTING, print Execution Engineer invocation
    resume         — Auto-detect current phase and print the next invocation
    recover        — Generate a rich RESUME_PROMPT.md and copy it to clipboard
    loop           — Print the full loop diagram and current position
    config         — Show the active configuration loaded from research_config.json

.PARAMETER StudyName
    The study folder name under studies/

.PARAMETER StudyGoal
    The research goal (required for init)

.PARAMETER RunDir
    Task run directory. Defaults to task_runs/<study_name>
#>
param(
    [ValidateSet('quickstart','init','status','plan','execute','review','validate','report','resume','recover','loop','config')]
    [string]$Action = 'status',

    [Parameter(Mandatory=$true)]
    [string]$StudyName,

    [string]$StudyGoal = '',

    [string]$RunDir = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$WorkspaceRoot = Split-Path -Parent $PSScriptRoot
$StudyDir      = Join-Path (Join-Path $WorkspaceRoot "studies") $StudyName
$ConfigFile    = Join-Path $WorkspaceRoot "research_config.json"

if ([string]::IsNullOrWhiteSpace($RunDir)) {
    $RunDir = Join-Path (Join-Path $WorkspaceRoot "task_runs") $StudyName
} elseif (-not [System.IO.Path]::IsPathRooted($RunDir)) {
    $RunDir = Join-Path $WorkspaceRoot $RunDir
}

$StudyStateFile      = Join-Path $StudyDir "study_state.json"
$ScienceDirective    = Join-Path $RunDir "SCIENCE_DIRECTIVE.md"
$ExecutionSummary    = Join-Path $RunDir "EXECUTION_SUMMARY.md"
$TaskChecklist       = Join-Path $RunDir "TASK_CHECKLIST.md"
$ProgressLog         = Join-Path $RunDir "PROGRESS_LOG.md"
$BlockersFile        = Join-Path $RunDir "BLOCKERS.md"
$DoneFile            = Join-Path $RunDir "DONE.md"
$RecoveryPromptFile  = Join-Path $RunDir "RESUME_PROMPT.md"
$FollowupFile        = Join-Path $RunDir "FOLLOWUP_PROMPT.md"
$ReviewRequestFile   = Join-Path $RunDir "REVIEW_REQUEST.md"
$ReviewDirectiveFile = Join-Path $RunDir "REVIEW_DIRECTIVE.md"
$PolishDoneFile      = Join-Path $RunDir "POLISH_COMPLETE.md"

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

function Get-RelPath {
    param([string]$AbsPath)
    $root = [System.IO.Path]::GetFullPath($WorkspaceRoot).TrimEnd('\','/')
    $abs  = [System.IO.Path]::GetFullPath($AbsPath)
    if ($abs.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
        return $abs.Substring($root.Length).TrimStart('\','/').Replace('\','/')
    }
    return $abs.Replace('\','/')
}

function Get-Timestamp {
    return (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ')
}

function Read-Config {
    if (Test-Path -LiteralPath $ConfigFile) {
        try {
            return Get-Content -LiteralPath $ConfigFile -Raw | ConvertFrom-Json
        } catch {
            Write-Warning "Could not parse research_config.json: $_. Using defaults."
        }
    }
    # Return built-in defaults when no config file exists
    return [PSCustomObject]@{
        models  = [PSCustomObject]@{
            planning      = [PSCustomObject]@{ model_id = "o3" }
            execution     = [PSCustomObject]@{ model_id = "claude-opus-4-6" }
            documentation = [PSCustomObject]@{ model_id = "claude-opus-4-6" }
        }
        loop    = [PSCustomObject]@{ max_iterations = 10; min_iterations = 1 }
        retry   = [PSCustomObject]@{ max_retries_per_phase = 3; auto_recover_on_failure = $true; blocked_phase_policy = "continue_with_partial" }
        report  = [PSCustomObject]@{ preserve_existing_report = $true; extension_mode = "append_iteration_section"; backup_before_write = $true }
        logging = [PSCustomObject]@{ generate_recovery_prompt = $true; recovery_prompt_file = "RESUME_PROMPT.md"; max_execution_summary_lines = 500 }
    }
}

function Get-MaxIterations {
    $cfg = Read-Config
    if ($cfg -and $cfg.loop -and $cfg.loop.max_iterations) { return [int]$cfg.loop.max_iterations }
    return 10
}

function Get-PlanningModel {
    $cfg = Read-Config
    if ($cfg -and $cfg.models -and @($cfg.models.PSObject.Properties.Name) -contains 'planning') {
        if ($cfg.models.planning -and $cfg.models.planning.model_id) {
            return $cfg.models.planning.model_id
        }
    }
    if ($cfg -and $cfg.models -and @($cfg.models.PSObject.Properties.Name) -contains 'review') {
        if ($cfg.models.review -and $cfg.models.review.model_id) {
            return $cfg.models.review.model_id
        }
    }
    return "o3"
}

function Get-ExecutionModel {
    $cfg = Read-Config
    if ($cfg -and $cfg.models -and @($cfg.models.PSObject.Properties.Name) -contains 'execution') {
        if ($cfg.models.execution -and $cfg.models.execution.model_id) {
            return $cfg.models.execution.model_id
        }
    }
    return "claude-opus-4-6"
}

function Read-StudyState {
    if (Test-Path -LiteralPath $StudyStateFile) {
        return Get-Content -LiteralPath $StudyStateFile -Raw | ConvertFrom-Json
    }
    return $null
}

function Write-StudyState {
    param([PSObject]$State)
    $State.updated_at = Get-Timestamp
    $json = $State | ConvertTo-Json -Depth 10
    $enc  = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($StudyStateFile, $json, $enc)
}

function Append-ProgressLog {
    param([string]$Entry)
    $ts   = Get-Timestamp
    $line = "`n## $ts - $Entry`n"
    Add-Content -LiteralPath $ProgressLog -Value $line -Encoding utf8
}

function Get-StateUpdatedAtUtc {
    param([PSObject]$State)

    if (-not $State -or -not $State.updated_at) {
        return $null
    }

    try {
        return [DateTime]::Parse(
            [string]$State.updated_at,
            [System.Globalization.CultureInfo]::InvariantCulture,
            [System.Globalization.DateTimeStyles]::AssumeUniversal -bor [System.Globalization.DateTimeStyles]::AdjustToUniversal
        )
    } catch {
        return $null
    }
}

function Test-IsExtensionPass {
    param([PSObject]$State)

    if (-not $State -or -not $State.objective) {
        return $false
    }

    return ([string]$State.objective) -match '\|\s*EXTENSION:'
}

function Get-IterationSignalFiles {
    return @(
        $ScienceDirective,
        $ExecutionSummary,
        $ReviewRequestFile,
        $ReviewDirectiveFile,
        $FollowupFile,
        $DoneFile,
        $PolishDoneFile
    )
}

function Move-StaleIterationSignals {
    param([PSObject]$State)

    $cutoff = Get-StateUpdatedAtUtc -State $State
    if (-not $cutoff) {
        return @()
    }

    $archiveRoot = Join-Path $RunDir "iteration_history"
    $archiveDir = $null
    $archived = @()

    foreach ($path in @(Get-IterationSignalFiles)) {
        if (-not (Test-Path -LiteralPath $path)) {
            continue
        }

        try {
            $item = Get-Item -LiteralPath $path -ErrorAction Stop
        } catch {
            continue
        }

        if ($item.LastWriteTimeUtc -ge $cutoff) {
            continue
        }

        if (-not $archiveDir) {
            if (-not (Test-Path -LiteralPath $archiveRoot)) {
                New-Item -ItemType Directory -Path $archiveRoot | Out-Null
            }
            $archiveDir = Join-Path $archiveRoot ("reset_" + $cutoff.ToString('yyyyMMddTHHmmssZ'))
            if (-not (Test-Path -LiteralPath $archiveDir)) {
                New-Item -ItemType Directory -Path $archiveDir | Out-Null
            }
        }

        Move-Item -LiteralPath $item.FullName -Destination (Join-Path $archiveDir $item.Name) -Force
        $archived += $item.Name
    }

    if (@($archived).Count -gt 0 -and (Test-Path -LiteralPath $ProgressLog)) {
        $ts = Get-Timestamp
        $lines = @("", "## $ts - Archived stale extension signal files")
        foreach ($name in @($archived)) {
            $lines += "- $name"
        }
        Add-Content -LiteralPath $ProgressLog -Value ($lines -join "`n") -Encoding utf8
    }

    return @($archived)
}

function Get-ActiveBlockers {
    if (-not (Test-Path -LiteralPath $BlockersFile)) { return @() }
    $blockers = @()
    $inActive = $false
    foreach ($ln in [System.IO.File]::ReadAllLines($BlockersFile)) {
        if ($ln -match '^## Active Blockers') { $inActive = $true; continue }
        if ($inActive -and $ln -match '^## ')  { break }
        if ($inActive -and $ln.Trim().StartsWith('-') -and $ln.Trim() -ne '- None.') {
            $blockers += $ln.Trim()
        }
    }
    return $blockers
}

function Get-OpenTasks {
    if (-not (Test-Path -LiteralPath $TaskChecklist)) { return @() }
    return @(Select-String -Path $TaskChecklist -Pattern '^- \[ \]' | ForEach-Object { $_.Line.Trim() })
}

function Get-CompletedTasks {
    if (-not (Test-Path -LiteralPath $TaskChecklist)) { return @() }
    return @(Select-String -Path $TaskChecklist -Pattern '^- \[x\]' | ForEach-Object { $_.Line.Trim() })
}

function Detect-CurrentPhase {
    $s = Read-StudyState
    if (-not $s) { return 'needs-init' }
    switch ($s.status) {
        'INITIALIZED'   { return 'needs-plan' }
        'PLANNING'      { return 'needs-plan' }
        'PLANNED'       { return 'needs-execute' }
        'IMPLEMENTING'  { return 'needs-execute' }
        'REVISION_REQUESTED' { return 'needs-execute' }
        'REVIEW_REQUESTED' { return 'needs-review' }
        'READY_FOR_REVIEW' { return 'needs-review' }
        'REVIEWING'     { return 'needs-review' }
        'VALIDATING'    { return 'needs-validate' }
        'REPORTING'     { return 'needs-report' }
        'APPROVED'      { return 'needs-polish' }
        'BLOCKED'       { return 'blocked' }
        'COMPLETE'      { return 'complete' }
        default         { return 'unknown' }
    }
}

function Write-PhaseHeader {
    param([string]$Phase, [string]$Description)
    Write-Output ""
    Write-Output ("=" * 64)
    Write-Output "  RESEARCH LOOP -- $Phase"
    Write-Output "  Study: $StudyName"
    Write-Output "  $Description"
    Write-Output ("=" * 64)
    Write-Output ""
}

# ---------------------------------------------------------------------------
# Recovery prompt generator
# ---------------------------------------------------------------------------
# Writes RESUME_PROMPT.md in the run directory and copies to clipboard.
# Includes the exact @agent invocation plus all context needed to continue.

function Write-RecoveryPrompt {
    param([PSObject]$State)
    if (-not $State) { return }

    $phase         = Detect-CurrentPhase
    $relStudy      = Get-RelPath $StudyDir
    $relRun        = Get-RelPath $RunDir
    $planModel     = Get-PlanningModel
    $execModel     = Get-ExecutionModel
    $maxIter       = Get-MaxIterations
    $curIter       = [int]$State.loop_iteration
    $ts            = Get-Timestamp
    $activeBlockers = @(Get-ActiveBlockers)
    $openTasks     = @(Get-OpenTasks)
    $cfg           = Read-Config
    $hasFollowup   = Test-Path -LiteralPath $FollowupFile

    # Determine which agent and model to invoke
    $agentInvoke = switch ($phase) {
        'needs-plan'     { "@science-director study=$relStudy run=$relRun phase=plan" }
        'needs-execute'  { "@execution-engineer study=$relStudy run=$relRun phase=implement" }
        'needs-review'   { "@science-director study=$relStudy run=$relRun phase=review" }
        'needs-validate' { "@execution-engineer study=$relStudy run=$relRun phase=validate" }
        'needs-report'   { "@execution-engineer study=$relStudy run=$relRun phase=report" }
        'needs-polish'   { "@execution-engineer study=$relStudy run=$relRun phase=polish" }
        'blocked'        { "@research-loop study=$relStudy resume   # Study is BLOCKED - review blockers first" }
        default          { "@research-loop study=$relStudy resume" }
    }

    $agentLabel = switch ($phase) {
        'needs-plan'     { "science-director  [model: $planModel]" }
        'needs-execute'  { "execution-engineer  [model: $execModel]" }
        'needs-review'   { "science-director  [model: $planModel]" }
        'needs-validate' { "execution-engineer  [model: $execModel]" }
        'needs-report'   { "execution-engineer  [model: $execModel]" }
        'needs-polish'   { "execution-engineer  [model: $execModel]" }
        default          { "research-loop  [model: $execModel]" }
    }

    $blockersText = if ($activeBlockers.Count -gt 0) {
        "## Active Blockers`n" + ($activeBlockers -join "`n") + "`n"
    } else {
        "## Active Blockers`n- None`n"
    }

    $openTasksText = if ($openTasks.Count -gt 0) {
        $listed = ($openTasks | Select-Object -First 10) -join "`n"
        "## Open Tasks (first 10 of $($openTasks.Count))`n$listed`n"
    } elseif ($phase -eq 'needs-execute' -and $hasFollowup) {
        "## Open Tasks`n- Review requested additional work. Read ``$relRun/FOLLOWUP_PROMPT.md`` for the ordered next actions.`n"
    } else {
        "## Open Tasks`n- None remaining`n"
    }

    $reportPreserveNote = if ($cfg -and $cfg.report -and $cfg.report.preserve_existing_report) {
        "Report preservation: ENABLED (append_iteration_section - do NOT overwrite report.tex)"
    } else {
        "Report preservation: DISABLED (overwrite allowed)"
    }

    $content = @"
# Recovery Prompt - $StudyName
Generated: $ts

---

## How to Use This File

1. Copy the **Agent Invocation** below and paste it into Copilot Chat.
2. The agent will read the state files listed under **Context Files** and
   continue from exactly where the loop left off.
3. Do NOT redo work that is already marked complete in TASK_CHECKLIST.md.

---

## Current Loop State

| Field            | Value                          |
|------------------|-------------------------------|
| Study            | $($State.study_name)          |
| Status           | $($State.status)              |
| Phase            | $phase                        |
| Loop iteration   | $curIter / $maxIter           |
| Objective        | $($State.objective)           |

---

## Agent to Invoke

**Agent:** $agentLabel

Paste this into Copilot Chat:

```
$agentInvoke
```

---

## Context Files (read these first)

The agent MUST read these files before doing anything:

1. ``$relStudy/study_state.json`` - machine-readable study state
2. ``$relRun/SCIENCE_DIRECTIVE.md`` - last planning directive
3. ``$relRun/EXECUTION_SUMMARY.md`` - last execution results (if exists)
4. ``$relRun/TASK_CHECKLIST.md`` - which tasks are done / open
5. ``$relRun/BLOCKERS.md`` - known blockers
6. ``$relStudy/IMPROVEMENTS.md`` - improvement log

---

## Recovery Instructions for the Agent

This is a **recovery invocation** for a study that was interrupted or needs continuation.

- Current phase: **$phase**
- Iteration: $curIter (max allowed: $maxIter)
- $reportPreserveNote

Do NOT restart from scratch. Do NOT redo completed tasks. Read the context files
above, orient yourself to the current state, and continue the study forward.

$blockersText
$openTasksText

---

## Key Results Achieved So Far

$(
    if ($State.key_results -and $State.key_results.PSObject.Properties.Count -gt 0) {
        ($State.key_results.PSObject.Properties | ForEach-Object { "- **$($_.Name):** $($_.Value)" }) -join "`n"
    } else {
        "- No key results recorded in study_state.json yet."
    }
)

## Recent Progress (last 10 entries)

$(
    $progFile = Join-Path $RunDir "PROGRESS_LOG.md"
    if (Test-Path -LiteralPath $progFile) {
        $lines = @(Get-Content -LiteralPath $progFile -Encoding UTF8 -ErrorAction SilentlyContinue)
        if ($lines.Count -gt 10) {
            ($lines | Select-Object -Last 10) -join "`n"
        } elseif ($lines.Count -gt 0) {
            $lines -join "`n"
        } else {
            "- Progress log is empty."
        }
    } else {
        "- No PROGRESS_LOG.md found."
    }
)

---

## Loop Configuration (from research_config.json)

- Planning model:   $planModel
- Execution model:  $execModel
- Max iterations:   $maxIter
- Blocked policy:   $(if ($cfg -and $cfg.retry) { $cfg.retry.blocked_phase_policy } else { "continue_with_partial" })
- Report mode:      $(if ($cfg -and $cfg.report) { $cfg.report.extension_mode } else { "append_iteration_section" })

---
*Generated by tools/research_loop.ps1 -Action recover -StudyName $StudyName*
"@

    $enc = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($RecoveryPromptFile, $content, $enc)

    # Copy to clipboard if available
    try {
        $agentInvoke | Set-Clipboard
        Write-Output "  Clipboard: agent invocation copied."
    } catch {
        # Clipboard not available in all environments - silently skip
    }
}

# ---------------------------------------------------------------------------
# Actions
# ---------------------------------------------------------------------------

switch ($Action) {

    # -----------------------------------------------------------------------
    'quickstart' {
        Write-PhaseHeader "QUICKSTART" "One command to start or extend a study"

        $relStudy  = Get-RelPath $StudyDir
        $relRun    = Get-RelPath $RunDir
        $execModel = Get-ExecutionModel
        $isNew     = -not (Test-Path -LiteralPath $StudyStateFile)
        $isExtend  = $false

        # ---- Bootstrap if the study does not exist yet ----
        if ($isNew) {
            if ([string]::IsNullOrWhiteSpace($StudyGoal)) {
                throw 'StudyGoal is required when starting a new study. Use -StudyGoal "your research question"'
            }
            Write-Output "Study does not exist. Initializing..."
            Write-Output ""

            foreach ($d in @('scripts','data','figures','report')) {
                New-Item -ItemType Directory -Path (Join-Path $StudyDir $d) -Force | Out-Null
            }
            New-Item -ItemType Directory -Path $RunDir -Force | Out-Null

            $ts = Get-Timestamp
            $initState = [PSCustomObject]@{
                study_name                = $StudyName
                study_path                = "studies/$StudyName"
                status                    = 'INITIALIZED'
                problem_class             = @()
                created_at                = $ts
                updated_at                = $ts
                loop_iteration            = 0
                phase_retry_counts        = [PSCustomObject]@{}
                objective                 = $StudyGoal
                hypotheses                = @()
                assumptions               = @()
                success_criteria          = [PSCustomObject]@{}
                completed_tasks           = @()
                failed_tasks              = @()
                pending_tasks             = @()
                blocked_tasks             = @()
                key_results               = [PSCustomObject]@{}
                latest_figures            = @()
                blockers                  = @()
                compute_notes             = [PSCustomObject]@{}
                science_directive_version = 0
                file_manifest             = [PSCustomObject]@{
                    scripts = @()
                    data    = @()
                    figures = @()
                    report  = ''
                }
            }
            Write-StudyState $initState

            $enc = New-Object System.Text.UTF8Encoding($false)
            if (-not (Test-Path (Join-Path $StudyDir "README.md"))) {
                $goalForReadme = $StudyGoal
                if ($goalForReadme.Length -gt 500) {
                    $goalForReadme = $goalForReadme.Substring(0, 500).TrimEnd() + "`n`n(See full prompt in the prompt file.)"
                }
                $readme = "# $StudyName`n`n## Problem Class`n<!-- OPT | REP | DES | ANA -->`n`n## Motivation`n$goalForReadme`n`n## Goals`n`n## Methods`n`n## Expected Outcomes`n`n## Known Limitations`n`n## Status`nACTIVE`n"
                [System.IO.File]::WriteAllText((Join-Path $StudyDir "README.md"), $readme, $enc)
            }
            if (-not (Test-Path (Join-Path $StudyDir "IMPROVEMENTS.md"))) {
                $improvements = "# Improvement Log: $StudyName`n`n> Written for future agents.`n`n## Critical Gaps (P1)`n`n## Recommended Improvements (P2)`n`n## Open Questions`n`n## What Was Tried and Did Not Work`n`n## Resolved`n"
                [System.IO.File]::WriteAllText((Join-Path $StudyDir "IMPROVEMENTS.md"), $improvements, $enc)
            }
            if (-not (Test-Path $ProgressLog)) {
                $goalSummary = $StudyGoal
                if ($goalSummary.Length -gt 200) {
                    $goalSummary = $goalSummary.Substring(0, 200).TrimEnd() + "... (truncated, see prompt file)"
                }
                $plog = "# Progress Log`n`n## $ts - Study initialized via quickstart`n- Objective: $goalSummary`n- Study path: $relStudy`n- Run path:   $relRun`n"
                [System.IO.File]::WriteAllText($ProgressLog, $plog, $enc)
            }
            if (-not (Test-Path $BlockersFile)) {
                [System.IO.File]::WriteAllText($BlockersFile, "# Blockers`n`n## Active Blockers`n- None.`n`n## Resolved Blockers`n- None.`n", $enc)
            }
            if (-not (Test-Path $TaskChecklist)) {
                $cl = "# Task Checklist`n`n## Bootstrap`n- [x] B0.1 Initialize study directory`n- [ ] B0.2 Science Director produces first SCIENCE_DIRECTIVE.md`n`n## Phase 1: Planning`n`n## Phase 2: Implementation`n`n## Phase 3: Validation`n`n## Phase 4: Reporting`n"
                [System.IO.File]::WriteAllText($TaskChecklist, $cl, $enc)
            }
            $state = $initState
            Write-Output "  Created: $relStudy"
            Write-Output "  Created: $relRun"

        } else {
            # ---- Study exists — determine if we are resuming or extending ----
            $state = Read-StudyState
            $archivedSignals = @()

            if (Test-IsExtensionPass -State $state) {
                $archivedSignals = @(Move-StaleIterationSignals -State $state)
                if (@($archivedSignals).Count -gt 0) {
                    Write-Output "Recovered extension state by archiving stale stage signals..."
                    Write-Output "  Archived: $($archivedSignals -join ', ')"
                }
            }

            $phase = Detect-CurrentPhase

            if ($phase -eq 'complete') {
                $isExtend = $true
                # Mark as INITIALIZED again so the full loop re-runs on new material
                $state.status = 'INITIALIZED'
                if (-not [string]::IsNullOrWhiteSpace($StudyGoal)) {
                    # Append the new goal to the objective
                    $state.objective = $state.objective + " | EXTENSION: " + $StudyGoal
                }
                Write-StudyState $state
                $archivedSignals = @(Move-StaleIterationSignals -State $state)
                Write-Output "Study is COMPLETE. Preparing extension pass..."
                Write-Output "  Prior iteration count: $([int]$state.loop_iteration)"
                Write-Output "  Report will be EXTENDED (preserve_existing_report = true)"
                if (@($archivedSignals).Count -gt 0) {
                    Write-Output "  Archived stale stage signals: $($archivedSignals -join ', ')"
                }
            } else {
                Write-Output "Study exists and is in progress. Resuming from phase: $phase"
                Write-RecoveryPrompt $state
            }
        }

        # ---- Build the single Copilot command ----
        $state    = Read-StudyState
        $goalArg  = if (-not [string]::IsNullOrWhiteSpace($StudyGoal)) { " goal=`"$StudyGoal`"" } else { "" }
        $modeTag  = if ($isExtend) { " # extension pass" } elseif ($isNew) { " # new study" } else { " # resume" }

        $copilotCmd = if ($isNew -or $isExtend) {
            "@research-loop study=$relStudy$goalArg$modeTag"
        } else {
            "@research-loop study=$relStudy resume$modeTag"
        }

        # Write it to RESUME_PROMPT.md and clipboard
        Write-RecoveryPrompt $state

        Write-Output ""
        Write-Output ("=" * 64)
        Write-Output "  PASTE THIS INTO COPILOT CHAT:"
        Write-Output ("=" * 64)
        Write-Output ""
        Write-Output "  $copilotCmd"
        Write-Output ""
        Write-Output ("=" * 64)
        Write-Output ""
        Write-Output "Model to select:  $execModel  (or your preferred model)"
        Write-Output "Context file:     $relRun/RESUME_PROMPT.md"
        Write-Output ""

        try {
            $copilotCmd | Set-Clipboard
            Write-Output "Command is in your clipboard. Ctrl+V into Copilot Chat."
        } catch {
            Write-Output "Copy the command above manually into Copilot Chat."
        }

        Write-Output ""
        if ($isNew) {
            Write-Output "What the agent will do:"
            Write-Output "  1. Read research_config.json and the study README"
            Write-Output "  2. Look up the cqed_sim API (cqed-sim-lookup skill)"
            Write-Output "  3. Design experiments (Science Director hat)"
            Write-Output "  4. Write and run all simulation code (Execution Engineer hat)"
            Write-Output "  5. Review results and iterate until convergence"
            Write-Output "  6. Validate and write the final LaTeX report"
        } elseif ($isExtend) {
            Write-Output "What the agent will do:"
            Write-Output "  1. Read all prior study state and the existing report"
            Write-Output "  2. Design extension experiments based on IMPROVEMENTS.md"
            Write-Output "  3. Implement, run, and validate new results"
            Write-Output "  4. EXTEND report.tex with a new section (never overwrite)"
        } else {
            Write-Output "What the agent will do:"
            Write-Output "  1. Read RESUME_PROMPT.md to restore context"
            Write-Output "  2. Continue from the last incomplete phase"
            Write-Output "  3. Skip all tasks already marked done"
        }
    }

    # -----------------------------------------------------------------------
    'init' {
        Write-PhaseHeader "INITIALIZE" "Creating study structure and task run"

        if ([string]::IsNullOrWhiteSpace($StudyGoal)) {
            throw 'StudyGoal is required for init. Use -StudyGoal "your research question"'
        }

        foreach ($d in @('scripts','data','figures','report')) {
            New-Item -ItemType Directory -Path (Join-Path $StudyDir $d) -Force | Out-Null
        }
        New-Item -ItemType Directory -Path $RunDir -Force | Out-Null

        $ts = Get-Timestamp
        $initState = [PSCustomObject]@{
            study_name                = $StudyName
            study_path                = "studies/$StudyName"
            status                    = 'INITIALIZED'
            problem_class             = @()
            created_at                = $ts
            updated_at                = $ts
            loop_iteration            = 0
            phase_retry_counts        = [PSCustomObject]@{}
            objective                 = $StudyGoal
            hypotheses                = @()
            assumptions               = @()
            success_criteria          = [PSCustomObject]@{}
            completed_tasks           = @()
            failed_tasks              = @()
            pending_tasks             = @()
            blocked_tasks             = @()
            key_results               = [PSCustomObject]@{}
            latest_figures            = @()
            blockers                  = @()
            compute_notes             = [PSCustomObject]@{}
            science_directive_version = 0
            file_manifest             = [PSCustomObject]@{
                scripts = @()
                data    = @()
                figures = @()
                report  = ''
            }
        }
        Write-StudyState $initState

        $enc      = New-Object System.Text.UTF8Encoding($false)
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir

        if (-not (Test-Path (Join-Path $StudyDir "README.md"))) {
            $readme = "# $StudyName`n`n## Problem Class`n<!-- OPT | REP | DES | ANA -->`n`n## Motivation`n$StudyGoal`n`n## Goals`n<!-- Numbered, concrete, falsifiable goals. -->`n`n## Methods`n<!-- Which cqed_sim modules/functions will be used. -->`n`n## Expected Outcomes`n<!-- What success looks like. -->`n`n## Known Limitations`n<!-- Updated throughout the study. -->`n`n## Status`nACTIVE`n"
            [System.IO.File]::WriteAllText((Join-Path $StudyDir "README.md"), $readme, $enc)
        }

        if (-not (Test-Path (Join-Path $StudyDir "IMPROVEMENTS.md"))) {
            $improvements = "# Improvement Log: $StudyName`n`n> This file is written for future agents. Be specific, honest, and actionable.`n`n## Critical Gaps (P1)`n`n## Recommended Improvements (P2)`n`n## Nice-to-Haves (P3)`n`n## Open Questions`n`n## What Was Tried and Did Not Work`n`n## Compute & Resource Notes`n`n## Resolved`n"
            [System.IO.File]::WriteAllText((Join-Path $StudyDir "IMPROVEMENTS.md"), $improvements, $enc)
        }

        if (-not (Test-Path $ProgressLog)) {
            $plog = "# Progress Log`n`n## $ts - Study initialized`n- Objective: $StudyGoal`n- Study path: $relStudy`n- Run path:   $relRun`n- Next: Science Director planning phase`n"
            [System.IO.File]::WriteAllText($ProgressLog, $plog, $enc)
        }

        if (-not (Test-Path $BlockersFile)) {
            $blk = "# Blockers`n`n## Active Blockers`n- None.`n`n## Resolved Blockers`n- None.`n"
            [System.IO.File]::WriteAllText($BlockersFile, $blk, $enc)
        }

        if (-not (Test-Path $TaskChecklist)) {
            $cl = "# Task Checklist`n`n## Status Summary`n- Study: $relStudy`n- Run:   $relRun`n- Loop iteration: 0`n`n## Bootstrap`n- [x] B0.1 Initialize study directory and state files`n- [ ] B0.2 Science Director produces first SCIENCE_DIRECTIVE.md`n`n## Phase 1: Planning`n`n## Phase 2: Implementation`n`n## Phase 3: Validation`n`n## Phase 4: Reporting`n"
            [System.IO.File]::WriteAllText($TaskChecklist, $cl, $enc)
        }

        # Write initial recovery prompt
        Write-RecoveryPrompt $initState

        Write-Output "Study initialized:"
        Write-Output "  Study directory:  $relStudy"
        Write-Output "  Run directory:    $relRun"
        Write-Output "  State file:       $(Get-RelPath $StudyStateFile)"
        Write-Output "  Config file:      $(Get-RelPath $ConfigFile)"
        Write-Output ""
        Write-Output "Planning model:   $(Get-PlanningModel)"
        Write-Output "Execution model:  $(Get-ExecutionModel)"
        Write-Output "Max iterations:   $(Get-MaxIterations)"
        Write-Output ""
        Write-Output "Next step - invoke Science Director in Copilot Chat:"
        Write-Output "  @science-director study=$relStudy run=$relRun phase=plan"
    }

    # -----------------------------------------------------------------------
    'status' {
        Write-PhaseHeader "STATUS" "Current state of the research loop"

        $state = Read-StudyState
        if (-not $state) {
            Write-Output "No study_state.json found. Run init first."
            return
        }

        $openTasks      = @(Get-OpenTasks)
        $completedTasks = @(Get-CompletedTasks)
        $activeBlockers = @(Get-ActiveBlockers)
        $phase          = Detect-CurrentPhase
        $maxIter        = Get-MaxIterations
        $curIter        = [int]$state.loop_iteration
        $relStudy       = Get-RelPath $StudyDir
        $relRun         = Get-RelPath $RunDir

        Write-Output "Study:            $($state.study_name)"
        Write-Output "Objective:        $($state.objective)"
        Write-Output "Status:           $($state.status)"
        Write-Output "Loop iteration:   $curIter / $maxIter"
        Write-Output "Current phase:    $phase"
        $probClass = if ($state.problem_class) { $state.problem_class -join ', ' } else { '(not set)' }
        Write-Output "Problem class:    $probClass"
        Write-Output ""
        Write-Output "Tasks completed:  $(@($completedTasks).Count)"
        Write-Output "Tasks open:       $(@($openTasks).Count)"
        Write-Output "Active blockers:  $(@($activeBlockers).Count)"
        Write-Output ""
        Write-Output "Models:"
        Write-Output "  Planning (Science Director):   $(Get-PlanningModel)"
        Write-Output "  Execution (Execution Engineer): $(Get-ExecutionModel)"
        Write-Output ""

        $hasDone      = Test-Path -LiteralPath $DoneFile
        $hasDirective = Test-Path -LiteralPath $ScienceDirective
        $hasSummary   = Test-Path -LiteralPath $ExecutionSummary
        $hasRecovery  = Test-Path -LiteralPath $RecoveryPromptFile

        Write-Output "Files:"
        Write-Output "  study_state.json:      EXISTS"
        if ($hasDirective) { Write-Output "  SCIENCE_DIRECTIVE.md:  EXISTS" } else { Write-Output "  SCIENCE_DIRECTIVE.md:  MISSING" }
        if ($hasSummary)   { Write-Output "  EXECUTION_SUMMARY.md:  EXISTS" } else { Write-Output "  EXECUTION_SUMMARY.md:  MISSING" }
        if ($hasRecovery)  { Write-Output "  RESUME_PROMPT.md:      EXISTS" } else { Write-Output "  RESUME_PROMPT.md:      (not yet generated)" }
        if ($hasDone)      { Write-Output "  DONE.md:               EXISTS" } else { Write-Output "  DONE.md:               MISSING" }
        Write-Output ""

        if (@($activeBlockers).Count -gt 0) {
            Write-Output "Active Blockers:"
            $activeBlockers | ForEach-Object { Write-Output "  $_" }
            Write-Output ""
        }

        if (@($openTasks).Count -gt 0) {
            Write-Output "Next open tasks:"
            $openTasks | Select-Object -First 5 | ForEach-Object { Write-Output "  $_" }
            Write-Output ""
        }

        Write-Output "=== Suggested Next Action ==="
        switch ($phase) {
            'needs-init'     { Write-Output "  .\research_loop.ps1 -Action init -StudyName $StudyName -StudyGoal '<goal>'" }
            'needs-plan'     { Write-Output "  @science-director study=$relStudy run=$relRun phase=plan" }
            'needs-execute'  { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=implement" }
            'needs-review'   { Write-Output "  @science-director study=$relStudy run=$relRun phase=review" }
            'needs-validate' { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=validate" }
            'needs-report'   { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=report" }
            'blocked'        { Write-Output "  Study is BLOCKED. Resolve blockers, then run: .\research_loop.ps1 -Action recover -StudyName $StudyName" }
            'complete'       { Write-Output "  Study is COMPLETE. Report: $relStudy/report/report.pdf" }
            default          { Write-Output "  Unknown state. Check study_state.json manually." }
        }
        Write-Output ""
        Write-Output "To generate a recovery prompt: .\research_loop.ps1 -Action recover -StudyName $StudyName"
    }

    # -----------------------------------------------------------------------
    'plan' {
        Write-PhaseHeader "PLAN" "Science Director planning phase"

        $state   = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }

        $curIter = [int]$state.loop_iteration
        $maxIter = Get-MaxIterations

        if ($curIter -ge $maxIter) {
            Write-Output "WARNING: Iteration limit reached ($curIter / $maxIter)."
            Write-Output "Skipping further planning. Transitioning to VALIDATING."
            Write-Output "To allow more iterations, increase max_iterations in research_config.json."
            $state.status = 'VALIDATING'
            Write-StudyState $state
            Write-RecoveryPrompt $state
            return
        }

        $state.status = 'PLANNING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-RecoveryPrompt $state

        Write-Output "State: PLANNING  (iteration $curIter / $maxIter)"
        Write-Output "Planning model: $(Get-PlanningModel)"
        Write-Output ""
        Write-Output "Invoke the Science Director in Copilot Chat:"
        Write-Output "  @science-director study=$relStudy run=$relRun phase=plan"
    }

    # -----------------------------------------------------------------------
    'execute' {
        Write-PhaseHeader "EXECUTE" "Execution Engineer implementation phase"

        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        if (-not (Test-Path -LiteralPath $ScienceDirective)) {
            throw "No SCIENCE_DIRECTIVE.md found. Run the plan phase first."
        }

        $state.status = 'IMPLEMENTING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-RecoveryPrompt $state

        Write-Output "State: IMPLEMENTING  (iteration $([int]$state.loop_iteration) / $(Get-MaxIterations))"
        Write-Output "Execution model: $(Get-ExecutionModel)"
        Write-Output ""
        Write-Output "Invoke the Execution Engineer in Copilot Chat:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=implement"
    }

    # -----------------------------------------------------------------------
    'review' {
        Write-PhaseHeader "REVIEW" "Science Director review phase"

        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        if (-not (Test-Path -LiteralPath $ExecutionSummary)) {
            throw "No EXECUTION_SUMMARY.md found. Run the execute phase first."
        }

        $curIter = [int]$state.loop_iteration + 1
        $maxIter = Get-MaxIterations

        $state.status         = 'REVIEWING'
        $state.loop_iteration = $curIter
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-RecoveryPrompt $state

        if ($curIter -ge $maxIter) {
            Write-Output "NOTE: This is iteration $curIter / $maxIter. If Science Director decides CONTINUE,"
            Write-Output "      the loop will exceed max_iterations. It should prefer VALIDATE or STOP."
            Write-Output ""
        }

        Write-Output "State: REVIEWING  (iteration $curIter / $maxIter)"
        Write-Output "Planning model: $(Get-PlanningModel)"
        Write-Output ""
        Write-Output "Invoke the Science Director in Copilot Chat:"
        Write-Output "  @science-director study=$relStudy run=$relRun phase=review"
    }

    # -----------------------------------------------------------------------
    'validate' {
        Write-PhaseHeader "VALIDATE" "Execution Engineer validation phase"

        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }

        $state.status = 'VALIDATING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-RecoveryPrompt $state

        Write-Output "State: VALIDATING"
        Write-Output "Execution model: $(Get-ExecutionModel)"
        Write-Output ""
        Write-Output "Invoke the Execution Engineer in Copilot Chat:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=validate"
    }

    # -----------------------------------------------------------------------
    'report' {
        Write-PhaseHeader "REPORT" "Execution Engineer report phase"

        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }

        $cfg              = Read-Config
        $preserveReport   = if ($cfg -and $cfg.report) { [bool]$cfg.report.preserve_existing_report } else { $true }
        $extensionMode    = if ($cfg -and $cfg.report -and $cfg.report.extension_mode) { $cfg.report.extension_mode } else { "append_iteration_section" }

        $state.status = 'REPORTING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-RecoveryPrompt $state

        Write-Output "State: REPORTING"
        Write-Output "Execution model:   $(Get-ExecutionModel)"
        Write-Output "Report mode:       $extensionMode"
        Write-Output "Preserve existing: $preserveReport"
        Write-Output ""
        Write-Output "Invoke the Execution Engineer in Copilot Chat:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=report"

        if ($preserveReport) {
            $reportTex = Join-Path $StudyDir "report/report.tex"
            if (Test-Path -LiteralPath $reportTex) {
                Write-Output ""
                Write-Output "IMPORTANT: report.tex already exists."
                Write-Output "  Mode: $extensionMode - the agent must READ the existing report first and EXTEND it."
                Write-Output "  The agent must NOT overwrite prior content."
            }
        }
    }

    # -----------------------------------------------------------------------
    'resume' {
        Write-PhaseHeader "RESUME" "Detecting current phase and printing next action"

        $state    = Read-StudyState
        $phase    = Detect-CurrentPhase
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir

        Write-Output "Detected phase: $phase"
        Write-Output ""

        switch ($phase) {
            'needs-init'     { Write-Output "  .\research_loop.ps1 -Action init -StudyName $StudyName -StudyGoal '<goal>'" }
            'needs-plan'     { Write-Output "  @science-director study=$relStudy run=$relRun phase=plan" }
            'needs-execute'  { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=implement" }
            'needs-review'   { Write-Output "  @science-director study=$relStudy run=$relRun phase=review" }
            'needs-validate' { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=validate" }
            'needs-report'   { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=report" }
            'needs-polish'   { Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=polish" }
            'blocked'        {
                $bl = @(Get-ActiveBlockers)
                Write-Output "Study is BLOCKED. Active blockers:"
                $bl | ForEach-Object { Write-Output "  $_" }
                Write-Output ""
                Write-Output "Resolve blockers, then run:"
                Write-Output "  .\research_loop.ps1 -Action recover -StudyName $StudyName"
            }
            'complete' { Write-Output "Study is COMPLETE. Report: $relStudy/report/report.pdf" }
        }

        # Always regenerate recovery prompt on resume
        if ($state) {
            Write-RecoveryPrompt $state
            $relRecovery = Get-RelPath $RecoveryPromptFile
            Write-Output ""
            Write-Output "Recovery prompt updated: $relRecovery"
            Write-Output "Agent invocation is in your clipboard."
        }
    }

    # -----------------------------------------------------------------------
    'recover' {
        Write-PhaseHeader "RECOVER" "Generating rich recovery prompt for interrupted study"

        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }

        $phase    = Detect-CurrentPhase
        $relRun   = Get-RelPath $RunDir
        $relStudy = Get-RelPath $StudyDir

        Write-RecoveryPrompt $state

        Write-Output "Recovery prompt: $relRun/RESUME_PROMPT.md"
        Write-Output "Phase:           $phase"
        Write-Output "Agent:           $(switch ($phase) {
            'needs-plan'     { "science-director  [$(Get-PlanningModel)]" }
            'needs-execute'  { "execution-engineer  [$(Get-ExecutionModel)]" }
            'needs-review'   { "science-director  [$(Get-PlanningModel)]" }
            'needs-validate' { "execution-engineer  [$(Get-ExecutionModel)]" }
            'needs-report'   { "execution-engineer  [$(Get-ExecutionModel)]" }
            'needs-polish'   { "execution-engineer  [$(Get-ExecutionModel)]" }
            default { "research-loop" }
        })"
        Write-Output ""
        Write-Output "Next:"
        Write-Output "  1. Open $relRun/RESUME_PROMPT.md"
        Write-Output "  2. Copy the agent invocation and paste it into Copilot Chat."
        Write-Output "  3. The agent invocation is also in your clipboard."
        Write-Output ""
        Write-Output "Or use the VS Code task: 'Research: Generate Recovery Prompt'"
    }

    # -----------------------------------------------------------------------
    'loop' {
        Write-PhaseHeader "LOOP" "Full research loop diagram"

        $state    = Read-StudyState
        $phase    = Detect-CurrentPhase
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        $maxIter  = Get-MaxIterations
        $curIter  = if ($state) { [int]$state.loop_iteration } else { 0 }

        Write-Output "Loop diagram:"
        Write-Output "  1. PLAN    -> @science-director   phase=plan    [model: $(Get-PlanningModel)]"
        Write-Output "  2. EXECUTE -> @execution-engineer phase=implement [model: $(Get-ExecutionModel)]"
        Write-Output "  3. REVIEW  -> @science-director   phase=review  [model: $(Get-PlanningModel)]"
        Write-Output "  4.           Decision: CONTINUE/REVISE -> back to 2"
        Write-Output "  5. VALIDATE -> @execution-engineer phase=validate [model: $(Get-ExecutionModel)]"
        Write-Output "  6. REPORT   -> @execution-engineer phase=report   [model: $(Get-ExecutionModel)]"
        Write-Output ""
        Write-Output "Current position: $phase  (iteration $curIter / $maxIter)"
        Write-Output ""
        switch ($phase) {
            'needs-plan'    { Write-Output "Next: @science-director study=$relStudy run=$relRun phase=plan" }
            'needs-execute' { Write-Output "Next: @execution-engineer study=$relStudy run=$relRun phase=implement" }
            'needs-review'  { Write-Output "Next: @science-director study=$relStudy run=$relRun phase=review" }
            'needs-polish'  { Write-Output "Next: @execution-engineer study=$relStudy run=$relRun phase=polish" }
            default         { Write-Output "Next: .\research_loop.ps1 -Action resume -StudyName $StudyName" }
        }
    }

    # -----------------------------------------------------------------------
    'config' {
        Write-PhaseHeader "CONFIG" "Active research loop configuration"

        $cfg = Read-Config
        if (Test-Path -LiteralPath $ConfigFile) {
            Write-Output "Config file: $(Get-RelPath $ConfigFile)"
        } else {
            Write-Output "Config file: NOT FOUND (using built-in defaults)"
        }
        Write-Output ""
        Write-Output "Models:"
        Write-Output "  Planning (Science Director):    $($cfg.models.planning.model_id)"
        Write-Output "  Execution (Execution Engineer): $($cfg.models.execution.model_id)"
        Write-Output "  Documentation:                  $($cfg.models.documentation.model_id)"
        Write-Output ""
        Write-Output "Loop:"
        Write-Output "  Max iterations:        $($cfg.loop.max_iterations)"
        Write-Output "  Min iterations:        $($cfg.loop.min_iterations)"
        Write-Output ""
        Write-Output "Retry:"
        Write-Output "  Auto recover:          $($cfg.retry.auto_recover_on_failure)"
        Write-Output "  Max retries/phase:     $($cfg.retry.max_retries_per_phase)"
        Write-Output "  Blocked phase policy:  $($cfg.retry.blocked_phase_policy)"
        Write-Output ""
        Write-Output "Report:"
        Write-Output "  Preserve existing:     $($cfg.report.preserve_existing_report)"
        Write-Output "  Extension mode:        $($cfg.report.extension_mode)"
        Write-Output "  Backup before write:   $($cfg.report.backup_before_write)"
        Write-Output ""
        Write-Output "Logging:"
        Write-Output "  Generate recovery prompt: $($cfg.logging.generate_recovery_prompt)"
        Write-Output "  Recovery prompt file:     $($cfg.logging.recovery_prompt_file)"
        Write-Output "  Max exec summary lines:   $($cfg.logging.max_execution_summary_lines)"
    }
}
