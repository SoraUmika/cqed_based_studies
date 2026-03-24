<#
.SYNOPSIS
    Orchestrates the continuous two-model research loop for cQED studies.

.DESCRIPTION
    Manages the lifecycle of a research study through a Science Director
    and Execution Engineer loop. Handles bootstrapping, state management,
    phase transitions, and recovery from interrupted runs.

.PARAMETER Action
    The action to perform: init, status, plan, execute, review, validate, report, resume, loop

.PARAMETER StudyName
    The study folder name under studies/

.PARAMETER StudyGoal
    The research goal (required for init)

.PARAMETER RunDir
    Task run directory. Defaults to task_runs/<study_name>
#>
param(
    [ValidateSet('init','status','plan','execute','review','validate','report','resume','loop')]
    [string]$Action = 'status',

    [Parameter(Mandatory=$true)]
    [string]$StudyName,

    [string]$StudyGoal = '',

    [string]$RunDir = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$WorkspaceRoot = Split-Path -Parent $PSScriptRoot
$StudyDir = Join-Path (Join-Path $WorkspaceRoot "studies") $StudyName

if ([string]::IsNullOrWhiteSpace($RunDir)) {
    $RunDir = Join-Path (Join-Path $WorkspaceRoot "task_runs") $StudyName
} elseif (-not [System.IO.Path]::IsPathRooted($RunDir)) {
    $RunDir = Join-Path $WorkspaceRoot $RunDir
}

$StudyStateFile = Join-Path $StudyDir "study_state.json"
$ScienceDirective = Join-Path $RunDir "SCIENCE_DIRECTIVE.md"
$ExecutionSummary = Join-Path $RunDir "EXECUTION_SUMMARY.md"
$TaskChecklist = Join-Path $RunDir "TASK_CHECKLIST.md"
$ProgressLog = Join-Path $RunDir "PROGRESS_LOG.md"
$BlockersFile = Join-Path $RunDir "BLOCKERS.md"
$DoneFile = Join-Path $RunDir "DONE.md"

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
    $ts = Get-Timestamp
    $line = "`n## $ts - $Entry`n"
    Add-Content -LiteralPath $ProgressLog -Value $line -Encoding utf8
}

function Get-ActiveBlockers {
    if (-not (Test-Path -LiteralPath $BlockersFile)) { return @() }
    $blockers = @()
    $inActive = $false
    foreach ($ln in [System.IO.File]::ReadAllLines($BlockersFile)) {
        if ($ln -match '^## Active Blockers') { $inActive = $true; continue }
        if ($inActive -and $ln -match '^## ') { break }
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
        'REVIEWING'     { return 'needs-review' }
        'VALIDATING'    { return 'needs-validate' }
        'REPORTING'     { return 'needs-report' }
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

# -- Actions ---------------------------------------------------------------

switch ($Action) {

    'init' {
        Write-PhaseHeader "INITIALIZE" "Creating study structure and task run"

        if ([string]::IsNullOrWhiteSpace($StudyGoal)) {
            throw 'StudyGoal is required for init action. Use -StudyGoal "your research question"'
        }

        $dirs = @('scripts','data','figures','report')
        foreach ($d in $dirs) {
            New-Item -ItemType Directory -Path (Join-Path $StudyDir $d) -Force | Out-Null
        }
        New-Item -ItemType Directory -Path $RunDir -Force | Out-Null

        $timestamp = Get-Timestamp
        $initialState = [PSCustomObject]@{
            study_name                = $StudyName
            study_path                = "studies/$StudyName"
            status                    = 'INITIALIZED'
            problem_class             = @()
            created_at                = $timestamp
            updated_at                = $timestamp
            loop_iteration            = 0
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
        Write-StudyState $initialState

        $enc = New-Object System.Text.UTF8Encoding($false)
        $readmePath = Join-Path $StudyDir "README.md"
        if (-not (Test-Path -LiteralPath $readmePath)) {
            $readme = "# $StudyName`n`n## Problem Class`n<!-- OPT | REP | DES | ANA -->`n`n## Motivation`n$StudyGoal`n`n## Goals`n<!-- Numbered, concrete, falsifiable goals. -->`n`n## Methods`n<!-- Which cqed_sim modules/functions will be used. -->`n`n## Expected Outcomes`n<!-- What success looks like. -->`n`n## Known Limitations`n<!-- Updated throughout the study. -->`n`n## Status`nACTIVE`n"
            [System.IO.File]::WriteAllText($readmePath, $readme, $enc)
        }

        $improvePath = Join-Path $StudyDir "IMPROVEMENTS.md"
        if (-not (Test-Path -LiteralPath $improvePath)) {
            $improvements = "# Improvement Log: $StudyName`n`n> This file is written for future agents. Be specific, honest, and actionable.`n`n## Critical Gaps (P1)`n`n## Recommended Improvements (P2)`n`n## Nice-to-Haves (P3)`n`n## Open Questions`n`n## What Was Tried and Did Not Work`n`n## Compute & Resource Notes`n`n## Resolved`n"
            [System.IO.File]::WriteAllText($improvePath, $improvements, $enc)
        }

        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir

        if (-not (Test-Path (Join-Path $RunDir "PROGRESS_LOG.md"))) {
            $plog = "# Progress Log`n`n## $timestamp - Study initialized`n- Objective: $StudyGoal`n- Study path: $relStudy`n- Run path: $relRun`n- Next: Science Director planning phase`n"
            [System.IO.File]::WriteAllText($ProgressLog, $plog, $enc)
        }

        if (-not (Test-Path (Join-Path $RunDir "BLOCKERS.md"))) {
            $blk = "# Blockers`n`n## Active Blockers`n- None.`n`n## Resolved Blockers`n- None.`n"
            [System.IO.File]::WriteAllText($BlockersFile, $blk, $enc)
        }

        if (-not (Test-Path (Join-Path $RunDir "TASK_CHECKLIST.md"))) {
            $cl = "# Task Checklist`n`n## Status Summary`n- Study: $relStudy`n- Run: $relRun`n- Loop iteration: 0`n`n## Bootstrap`n- [x] B0.1 Initialize study directory and state files`n- [ ] B0.2 Science Director produces first SCIENCE_DIRECTIVE.md`n`n## Phase 1: Planning`n`n## Phase 2: Implementation`n`n## Phase 3: Validation`n`n## Phase 4: Reporting`n"
            [System.IO.File]::WriteAllText($TaskChecklist, $cl, $enc)
        }

        Write-Output "Study initialized:"
        Write-Output "  Study directory: $relStudy"
        Write-Output "  Run directory:   $relRun"
        Write-Output "  State file:      $(Get-RelPath $StudyStateFile)"
        Write-Output ""
        Write-Output "Next steps:"
        Write-Output "  1. @science-director study=$relStudy run=$relRun phase=plan"
        Write-Output "  2. Or: .\research_loop.ps1 -Action plan -StudyName $StudyName"
    }

    'status' {
        Write-PhaseHeader "STATUS" "Current state of the research loop"

        $state = Read-StudyState
        if (-not $state) {
            Write-Output "No study_state.json found. Run init first."
            return
        }

        $openTasks      = Get-OpenTasks
        $completedTasks = Get-CompletedTasks
        $activeBlockers = Get-ActiveBlockers
        $phase          = Detect-CurrentPhase

        Write-Output "Study:            $($state.study_name)"
        Write-Output "Objective:        $($state.objective)"
        Write-Output "Status:           $($state.status)"
        Write-Output "Loop iteration:   $($state.loop_iteration)"
        Write-Output "Current phase:    $phase"
        $probClass = if ($state.problem_class) { $state.problem_class -join ', ' } else { '(not set)' }
        Write-Output "Problem class:    $probClass"
        Write-Output ""
        Write-Output "Tasks completed:  $(@($completedTasks).Count)"
        Write-Output "Tasks open:       $(@($openTasks).Count)"
        Write-Output "Active blockers:  $(@($activeBlockers).Count)"
        Write-Output ""

        $hasDone      = Test-Path -LiteralPath $DoneFile
        $hasDirective = Test-Path -LiteralPath $ScienceDirective
        $hasSummary   = Test-Path -LiteralPath $ExecutionSummary

        Write-Output "Files:"
        Write-Output "  study_state.json:      EXISTS"
        if ($hasDirective) { Write-Output "  SCIENCE_DIRECTIVE.md:  EXISTS" } else { Write-Output "  SCIENCE_DIRECTIVE.md:  MISSING" }
        if ($hasSummary)   { Write-Output "  EXECUTION_SUMMARY.md:  EXISTS" } else { Write-Output "  EXECUTION_SUMMARY.md:  MISSING" }
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
            'needs-init'     { Write-Output "  Run: .\research_loop.ps1 -Action init -StudyName $StudyName -StudyGoal '<goal>'" }
            'needs-plan'     { Write-Output "  @science-director study=studies/$StudyName run=task_runs/$StudyName phase=plan" }
            'needs-execute'  { Write-Output "  @execution-engineer study=studies/$StudyName run=task_runs/$StudyName phase=implement" }
            'needs-review'   { Write-Output "  @science-director study=studies/$StudyName run=task_runs/$StudyName phase=review" }
            'needs-validate' { Write-Output "  @execution-engineer study=studies/$StudyName run=task_runs/$StudyName phase=validate" }
            'needs-report'   { Write-Output "  @execution-engineer study=studies/$StudyName run=task_runs/$StudyName phase=report" }
            'blocked'        { Write-Output "  Study is BLOCKED. Review blockers and resolve manually." }
            'complete'       { Write-Output "  Study is COMPLETE. Open report at studies/$StudyName/report/report.pdf" }
            default          { Write-Output "  Unknown state. Check study_state.json manually." }
        }
    }

    'plan' {
        Write-PhaseHeader "PLAN" "Science Director planning phase"
        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        $state.status = 'PLANNING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "State updated to PLANNING."
        Write-Output ""
        Write-Output "Invoke the Science Director agent with:"
        Write-Output "  @science-director study=$relStudy run=$relRun phase=plan"
    }

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
        Write-Output "State updated to IMPLEMENTING."
        Write-Output ""
        Write-Output "Invoke the Execution Engineer agent with:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=implement"
    }

    'review' {
        Write-PhaseHeader "REVIEW" "Science Director review phase"
        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        if (-not (Test-Path -LiteralPath $ExecutionSummary)) {
            throw "No EXECUTION_SUMMARY.md found. Run the execute phase first."
        }
        $state.status = 'REVIEWING'
        $state.loop_iteration = [int]$state.loop_iteration + 1
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "State updated to REVIEWING (iteration $($state.loop_iteration))."
        Write-Output ""
        Write-Output "Invoke the Science Director agent with:"
        Write-Output "  @science-director study=$relStudy run=$relRun phase=review"
    }

    'validate' {
        Write-PhaseHeader "VALIDATE" "Execution Engineer validation phase"
        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        $state.status = 'VALIDATING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "State updated to VALIDATING."
        Write-Output ""
        Write-Output "Invoke the Execution Engineer agent with:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=validate"
    }

    'report' {
        Write-PhaseHeader "REPORT" "Execution Engineer report phase"
        $state = Read-StudyState
        if (-not $state) { throw "No study state. Run init first." }
        $state.status = 'REPORTING'
        Write-StudyState $state
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "State updated to REPORTING."
        Write-Output ""
        Write-Output "Invoke the Execution Engineer agent with:"
        Write-Output "  @execution-engineer study=$relStudy run=$relRun phase=report"
    }

    'resume' {
        Write-PhaseHeader "RESUME" "Detecting current phase and continuing"
        $phase    = Detect-CurrentPhase
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "Detected phase: $phase"
        Write-Output ""
        switch ($phase) {
            'needs-init'     { Write-Output "  Run: .\research_loop.ps1 -Action init -StudyName $StudyName -StudyGoal '<goal>'" }
            'needs-plan'     { Write-Output "  Resuming -> @science-director study=$relStudy run=$relRun phase=plan" }
            'needs-execute'  { Write-Output "  Resuming -> @execution-engineer study=$relStudy run=$relRun phase=implement" }
            'needs-review'   { Write-Output "  Resuming -> @science-director study=$relStudy run=$relRun phase=review" }
            'needs-validate' { Write-Output "  Resuming -> @execution-engineer study=$relStudy run=$relRun phase=validate" }
            'needs-report'   { Write-Output "  Resuming -> @execution-engineer study=$relStudy run=$relRun phase=report" }
            'blocked'        {
                $bl = Get-ActiveBlockers
                Write-Output "Study is BLOCKED. Active blockers:"
                $bl | ForEach-Object { Write-Output "  $_" }
            }
            'complete' { Write-Output "  Study is COMPLETE. Report: studies/$StudyName/report/report.pdf" }
        }
    }

    'loop' {
        Write-PhaseHeader "LOOP" "Full research loop iteration"
        $phase    = Detect-CurrentPhase
        $relStudy = Get-RelPath $StudyDir
        $relRun   = Get-RelPath $RunDir
        Write-Output "Current phase: $phase"
        Write-Output ""
        Write-Output "The research loop runs:"
        Write-Output "  1. PLAN    -> @science-director  phase=plan"
        Write-Output "  2. EXECUTE -> @execution-engineer phase=implement"
        Write-Output "  3. REVIEW  -> @science-director  phase=review"
        Write-Output "  4. (loop to 2 if CONTINUE/REVISE)"
        Write-Output "  5. VALIDATE -> @execution-engineer phase=validate"
        Write-Output "  6. REPORT   -> @execution-engineer phase=report"
        Write-Output ""
        switch ($phase) {
            'needs-plan'    { Write-Output "Start: @science-director study=$relStudy run=$relRun phase=plan" }
            'needs-execute' { Write-Output "Start: @execution-engineer study=$relStudy run=$relRun phase=implement" }
            'needs-review'  { Write-Output "Start: @science-director study=$relStudy run=$relRun phase=review" }
            default         { Write-Output "Use: .\research_loop.ps1 -Action resume -StudyName $StudyName" }
        }
    }
}