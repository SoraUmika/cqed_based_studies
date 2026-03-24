[CmdletBinding()]
param(
    [ValidateSet('init', 'status')]
    [string]$Action = 'init',

    [string]$TaskFile = 'RESEARCH_PLAN.md',

    [string]$RunDir
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Get-WorkspaceRoot {
    return Split-Path -Parent $PSScriptRoot
}

function Resolve-WorkspacePath {
    param(
        [Parameter(Mandatory = $true)]
        [string]$PathValue,

        [Parameter(Mandatory = $true)]
        [string]$WorkspaceRoot
    )

    if ([System.IO.Path]::IsPathRooted($PathValue)) {
        return [System.IO.Path]::GetFullPath($PathValue)
    }

    return [System.IO.Path]::GetFullPath((Join-Path $WorkspaceRoot $PathValue))
}

function Get-WorkspaceRelativePath {
    param(
        [Parameter(Mandatory = $true)]
        [string]$AbsolutePath,

        [Parameter(Mandatory = $true)]
        [string]$WorkspaceRoot
    )

    $normalizedRoot = [System.IO.Path]::GetFullPath($WorkspaceRoot).TrimEnd([char[]]@('\', '/'))
    $normalizedPath = [System.IO.Path]::GetFullPath($AbsolutePath)

    if ($normalizedPath.StartsWith($normalizedRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
        return $normalizedPath.Substring($normalizedRoot.Length).TrimStart([char[]]@('\', '/')).Replace('\', '/')
    }

    return $normalizedPath.Replace('\', '/')
}

function Convert-ToRunSlug {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Value
    )

    $slug = [System.IO.Path]::GetFileNameWithoutExtension($Value).ToLowerInvariant()
    $slug = $slug -replace '[^a-z0-9]+', '_'
    $slug = $slug.Trim('_')

    if ([string]::IsNullOrWhiteSpace($slug)) {
        return 'task_run'
    }

    return $slug
}

function Write-FileIfMissing {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path,

        [Parameter(Mandatory = $true)]
        [string]$Content
    )

    if (Test-Path -LiteralPath $Path) {
        return $false
    }

    $utf8NoBom = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($Path, $Content, $utf8NoBom)
    return $true
}

function Get-ActiveBlockers {
    param(
        [Parameter(Mandatory = $true)]
        [string]$BlockersPath
    )

    if (-not (Test-Path -LiteralPath $BlockersPath)) {
        return @()
    }

    $activeBlockers = New-Object System.Collections.Generic.List[string]
    $inActiveSection = $false

    foreach ($line in [System.IO.File]::ReadAllLines($BlockersPath)) {
        if ($line -match '^## Active Blockers') {
            $inActiveSection = $true
            continue
        }

        if ($inActiveSection -and $line -match '^## ') {
            break
        }

        if ($inActiveSection) {
            $trimmed = $line.Trim()
            if ($trimmed.StartsWith('-') -and $trimmed -ne '- None.') {
                $activeBlockers.Add($trimmed)
            }
        }
    }

    return $activeBlockers.ToArray()
}

function Get-LastProgressCheckpoint {
    param(
        [Parameter(Mandatory = $true)]
        [string]$ProgressLogPath
    )

    if (-not (Test-Path -LiteralPath $ProgressLogPath)) {
        return $null
    }

    $entries = @(Select-String -Path $ProgressLogPath -Pattern '^## ' -AllMatches)
    if ($entries.Count -eq 0) {
        return $null
    }

    return $entries[$entries.Count - 1].Line.Trim()
}

$workspaceRoot = Get-WorkspaceRoot

if ([string]::IsNullOrWhiteSpace($RunDir)) {
    $RunDir = Join-Path 'task_runs' (Convert-ToRunSlug -Value $TaskFile)
}

$resolvedRunDir = Resolve-WorkspacePath -PathValue $RunDir -WorkspaceRoot $workspaceRoot
$relativeRunDir = Get-WorkspaceRelativePath -AbsolutePath $resolvedRunDir -WorkspaceRoot $workspaceRoot

if ($Action -eq 'init') {
    $resolvedTaskFile = Resolve-WorkspacePath -PathValue $TaskFile -WorkspaceRoot $workspaceRoot
    if (-not (Test-Path -LiteralPath $resolvedTaskFile -PathType Leaf)) {
        throw "Task file not found: $TaskFile"
    }

    $relativeTaskFile = Get-WorkspaceRelativePath -AbsolutePath $resolvedTaskFile -WorkspaceRoot $workspaceRoot
    $timestamp = (Get-Date).ToUniversalTime().ToString('yyyy-MM-ddTHH:mm:ssZ')

    New-Item -ItemType Directory -Path $resolvedRunDir -Force | Out-Null

    $executionPlanPath = Join-Path $resolvedRunDir 'EXECUTION_PLAN.md'
    $taskChecklistPath = Join-Path $resolvedRunDir 'TASK_CHECKLIST.md'
    $progressLogPath = Join-Path $resolvedRunDir 'PROGRESS_LOG.md'
    $blockersPath = Join-Path $resolvedRunDir 'BLOCKERS.md'

    $executionPlan = @"
# Execution Plan

## Run Metadata
- Source task: $relativeTaskFile
- Run directory: $relativeRunDir
- Status: NOT_PLANNED
- Last planner update: pending

## Task Summary
Pending planner synthesis from the source task document.

## Execution Rules
- Always re-read the source task and all run-state files before continuing work.
- Keep planning separate from implementation.
- Treat TASK_CHECKLIST.md as the execution source of truth.
- Append checkpoints to PROGRESS_LOG.md.
- Record unresolved issues in BLOCKERS.md.
- Write DONE.md only after final verification.

## Planned Phases
Planner will replace this section with ordered phases, dependencies, deliverables, and validation gates.

## Success Criteria
Planner will replace this section with explicit completion conditions.

## Continuation Protocol
1. Read the source task and every state file in this run directory.
2. Resume from the first unfinished checklist item that is not blocked.
3. Do not redo checked items unless a rollback is documented in PROGRESS_LOG.md.
4. Stop after reaching a meaningful checkpoint, encountering a blocker, or completing the task.
"@

    $taskChecklist = @"
# Task Checklist

## Status Summary
- Source task: $relativeTaskFile
- Run directory: $relativeRunDir
- Planner status: pending
- Implementation status: not started

## Bootstrap
- [x] R0.1 Create the run directory and bootstrap state files
- [ ] R0.2 Synthesize the execution plan from the source task
- [ ] R0.3 Complete all planned implementation checkpoints
- [ ] R0.4 Pass final verification and write DONE.md

## Planned Work
Planner will replace this section with phase-specific checklist items and stable task IDs.
"@

    $progressLog = @"
# Progress Log

## $timestamp - Run initialized
- Source task: $relativeTaskFile
- Run directory: $relativeRunDir
- Action: created bootstrap state files
- Next step: run /Autonomous Plan with task=$relativeTaskFile run=$relativeRunDir
"@

    $blockers = @"
# Blockers

## Active Blockers
- None.

## Resolved Blockers
- None.
"@

    $created = New-Object System.Collections.Generic.List[string]
    if (Write-FileIfMissing -Path $executionPlanPath -Content $executionPlan) { $created.Add('EXECUTION_PLAN.md') }
    if (Write-FileIfMissing -Path $taskChecklistPath -Content $taskChecklist) { $created.Add('TASK_CHECKLIST.md') }
    if (Write-FileIfMissing -Path $progressLogPath -Content $progressLog) { $created.Add('PROGRESS_LOG.md') }
    if (Write-FileIfMissing -Path $blockersPath -Content $blockers) { $created.Add('BLOCKERS.md') }

    Write-Output "Run directory: $relativeRunDir"
    Write-Output "Source task: $relativeTaskFile"
    if ($created.Count -gt 0) {
        Write-Output ('Created: ' + ($created -join ', '))
    }
    else {
        Write-Output 'Created: nothing new; state files already existed.'
    }
    Write-Output 'Next chat action: /Autonomous Plan'
    Write-Output "Suggested prompt text: task=$relativeTaskFile run=$relativeRunDir"
    return
}

if (-not (Test-Path -LiteralPath $resolvedRunDir -PathType Container)) {
    throw "Run directory not found: $RunDir"
}

$checklistPath = Join-Path $resolvedRunDir 'TASK_CHECKLIST.md'
$progressLogPath = Join-Path $resolvedRunDir 'PROGRESS_LOG.md'
$blockersPath = Join-Path $resolvedRunDir 'BLOCKERS.md'
$donePath = Join-Path $resolvedRunDir 'DONE.md'

$openTasks = @()
if (Test-Path -LiteralPath $checklistPath) {
    $openTasks = @(Select-String -Path $checklistPath -Pattern '^- \[ \]' | ForEach-Object { $_.Line.Trim() })
}

$activeBlockers = @(Get-ActiveBlockers -BlockersPath $blockersPath)
$lastCheckpoint = Get-LastProgressCheckpoint -ProgressLogPath $progressLogPath

Write-Output "Run directory: $relativeRunDir"
Write-Output ('DONE.md present: ' + $(if (Test-Path -LiteralPath $donePath) { 'yes' } else { 'no' }))
Write-Output ('Open checklist items: ' + $openTasks.Count)
Write-Output ('Active blockers: ' + $activeBlockers.Count)

if ($lastCheckpoint) {
    Write-Output ('Last checkpoint: ' + $lastCheckpoint)
}

if ($openTasks.Count -gt 0) {
    Write-Output 'Next unchecked items:'
    $openTasks | Select-Object -First 5 | ForEach-Object { Write-Output ('  ' + $_) }
}

if ($activeBlockers.Count -gt 0) {
    Write-Output 'Active blockers:'
    $activeBlockers | ForEach-Object { Write-Output ('  ' + $_) }
}