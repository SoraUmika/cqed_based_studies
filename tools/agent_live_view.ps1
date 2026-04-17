param(
    [Parameter(Mandatory=$true)]
    [string]$StudyName,

    [int]$PollSeconds = 3,

    [switch]$SnapshotOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$Root = Split-Path -Parent $PSScriptRoot
$RunDir = Join-Path (Join-Path $Root "task_runs") $StudyName
$StudyDir = Join-Path (Join-Path $Root "studies") $StudyName
$PolishDoneFile = Join-Path $RunDir "POLISH_COMPLETE.md"

if (-not (Test-Path -LiteralPath $RunDir)) {
    Write-Host "Run directory not found: $RunDir" -ForegroundColor Red
    exit 1
}

function Write-Banner {
    param([string]$Text, [string]$Color = "Cyan")
    $line = "=" * 70
    Write-Host ""
    Write-Host $line -ForegroundColor $Color
    Write-Host ("  " + $Text) -ForegroundColor $Color
    Write-Host $line -ForegroundColor $Color
    Write-Host ""
}

function Write-Step {
    param([string]$Text, [string]$Color = "White")
    $ts = (Get-Date).ToString("HH:mm:ss")
    Write-Host ("[" + $ts + "] " + $Text) -ForegroundColor $Color
}

function Get-TrackedFileEntries {
    return @(
        @{ Label = "PROGRESS_LOG";      Path = (Join-Path $RunDir "PROGRESS_LOG.md") },
        @{ Label = "TASK_CHECKLIST";    Path = (Join-Path $RunDir "TASK_CHECKLIST.md") },
        @{ Label = "BLOCKERS";          Path = (Join-Path $RunDir "BLOCKERS.md") },
        @{ Label = "SCIENCE_DIRECTIVE"; Path = (Join-Path $RunDir "SCIENCE_DIRECTIVE.md") },
        @{ Label = "EXECUTION_SUMMARY"; Path = (Join-Path $RunDir "EXECUTION_SUMMARY.md") },
        @{ Label = "REVIEW_REQUEST";    Path = (Join-Path $RunDir "REVIEW_REQUEST.md") },
        @{ Label = "REVIEW_DIRECTIVE";  Path = (Join-Path $RunDir "REVIEW_DIRECTIVE.md") },
        @{ Label = "FOLLOWUP_PROMPT";   Path = (Join-Path $RunDir "FOLLOWUP_PROMPT.md") },
        @{ Label = "DONE";              Path = (Join-Path $RunDir "DONE.md") },
        @{ Label = "POLISH_COMPLETE";   Path = (Join-Path $RunDir "POLISH_COMPLETE.md") },
        @{ Label = "STUDY_STATE";       Path = (Join-Path $StudyDir "study_state.json") }
    )
}

function Read-Lines {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        return @()
    }
    return @(Get-Content -LiteralPath $Path -Encoding UTF8)
}

function Show-LineBlock {
    param(
        [string]$Label,
        [string[]]$Lines,
        [string]$Color = "White"
    )

    Write-Host ""
    Write-Host (">>> " + $Label) -ForegroundColor Cyan
    foreach ($line in @($Lines)) {
        Write-Host $line -ForegroundColor $Color
    }
}

function Get-LastLines {
    param(
        [string[]]$Lines,
        [int]$MaxLines = 20
    )

    $count = @($Lines).Count
    if ($count -le $MaxLines) {
        return @($Lines)
    }
    return @($Lines | Select-Object -Last $MaxLines)
}

# ---------------------------------------------------------------------------
# Study directory file tracking (scripts, data, figures, report, artifacts)
# ---------------------------------------------------------------------------
$StudySubdirs = @('scripts', 'data', 'figures', 'report', 'artifacts')
$studyFileState = @{}   # path -> LastWriteTimeUtc

function Get-StudyFiles {
    $files = @()
    foreach ($sub in $StudySubdirs) {
        $subPath = Join-Path $StudyDir $sub
        if (Test-Path -LiteralPath $subPath) {
            $found = @(Get-ChildItem -LiteralPath $subPath -Recurse -File -ErrorAction SilentlyContinue)
            $files += $found
        }
    }
    # Also track top-level study files (README, IMPROVEMENTS, etc.)
    if (Test-Path -LiteralPath $StudyDir) {
        $topFiles = @(Get-ChildItem -LiteralPath $StudyDir -File -ErrorAction SilentlyContinue)
        $files += $topFiles
    }
    return $files
}

function Initialize-StudyFileState {
    $files = Get-StudyFiles
    foreach ($f in @($files)) {
        $studyFileState[$f.FullName] = $f.LastWriteTimeUtc
    }
    return @($files).Count
}

function Check-StudyFileChanges {
    $changes = @()
    $files = Get-StudyFiles
    $currentPaths = @{}

    foreach ($f in @($files)) {
        $currentPaths[$f.FullName] = $true
        $prior = $studyFileState[$f.FullName]
        if ($null -eq $prior) {
            # New file
            $relPath = $f.FullName
            if ($f.FullName.StartsWith($StudyDir)) {
                $relPath = $f.FullName.Substring($StudyDir.Length).TrimStart('\','/')
            }
            $changes += @{ Type = "NEW"; RelPath = $relPath; Size = $f.Length }
            $studyFileState[$f.FullName] = $f.LastWriteTimeUtc
        } elseif ($f.LastWriteTimeUtc -gt $prior) {
            # Modified file
            $relPath = $f.FullName
            if ($f.FullName.StartsWith($StudyDir)) {
                $relPath = $f.FullName.Substring($StudyDir.Length).TrimStart('\','/')
            }
            $changes += @{ Type = "MOD"; RelPath = $relPath; Size = $f.Length }
            $studyFileState[$f.FullName] = $f.LastWriteTimeUtc
        }
    }

    return $changes
}

function Format-FileSize {
    param([long]$Bytes)
    if ($Bytes -lt 1024) { return "$Bytes B" }
    if ($Bytes -lt 1048576) { return "{0:N1} KB" -f ($Bytes / 1024) }
    return "{0:N1} MB" -f ($Bytes / 1048576)
}

# ---------------------------------------------------------------------------
# Initialize state tracking
# ---------------------------------------------------------------------------
$trackedFiles = Get-TrackedFileEntries
$state = @{}
$lastHeartbeat = Get-Date
$HeartbeatSeconds = 300

Write-Banner "AGENT LIVE VIEW" "Green"
Write-Step ("Study: " + $StudyName) "White"
Write-Step ("Study dir: " + $StudyDir) "White"
Write-Step ("Run dir:   " + $RunDir) "White"
Write-Host ""
Write-Step "Monitors state files AND study directory for agent activity." "White"
Write-Step "For token-by-token live output, keep Copilot Chat open." "Yellow"
Write-Host ""

# Track state files (task_runs)
Write-Host "Tracked state files:" -ForegroundColor Green
foreach ($entry in $trackedFiles) {
    Write-Host ("  " + $entry.Label) -ForegroundColor DarkGray
    $lines = Read-Lines -Path $entry.Path
    $timestamp = [datetime]::MinValue
    if (Test-Path -LiteralPath $entry.Path) {
        $timestamp = (Get-Item -LiteralPath $entry.Path).LastWriteTimeUtc
    }
    $state[$entry.Path] = @{ Count = @($lines).Count; Timestamp = $timestamp }
}

# Track study directory files
$initialFileCount = Initialize-StudyFileState
Write-Host ""
Write-Host "Tracked study directories:" -ForegroundColor Green
foreach ($sub in $StudySubdirs) {
    $subPath = Join-Path $StudyDir $sub
    if (Test-Path -LiteralPath $subPath) {
        $subCount = @(Get-ChildItem -LiteralPath $subPath -Recurse -File -ErrorAction SilentlyContinue).Count
        Write-Host ("  " + $sub + "/ (" + $subCount + " files)") -ForegroundColor DarkGray
    } else {
        Write-Host ("  " + $sub + "/ (not yet created)") -ForegroundColor DarkGray
    }
}
Write-Host ("  Total tracked: " + $initialFileCount + " study files") -ForegroundColor DarkGray

# Show initial tails of state files
foreach ($entry in $trackedFiles) {
    $lines = Read-Lines -Path $entry.Path
    if (@($lines).Count -gt 0) {
        Show-LineBlock -Label ($entry.Label + " (initial tail)") -Lines (Get-LastLines -Lines $lines -MaxLines 12)
    }
}

if ($SnapshotOnly) {
    Write-Host "" 
    Write-Step "SnapshotOnly requested. Exiting." "Yellow"
    exit 0
}

Write-Host ""
Write-Step "Live monitor active. Ctrl+C to stop." "Green"
Write-Host ""

while ($true) {
    $anyChange = $false

    # ---- Check state files in task_runs ----
    foreach ($entry in $trackedFiles) {
        $path = $entry.Path
        $label = $entry.Label
        $prior = $state[$path]
        $lines = Read-Lines -Path $path
        $count = @($lines).Count
        $timestamp = [datetime]::MinValue
        if (Test-Path -LiteralPath $path) {
            $timestamp = (Get-Item -LiteralPath $path).LastWriteTimeUtc
        }

        $changed = ($timestamp -gt $prior.Timestamp) -or ($count -ne $prior.Count)
        if (-not $changed) {
            continue
        }

        $anyChange = $true

        if ($count -gt $prior.Count -and $prior.Count -gt 0) {
            $newLines = @($lines | Select-Object -Skip $prior.Count)
            if (@($newLines).Count -gt 0) {
                Show-LineBlock -Label ($label + " (new lines)") -Lines $newLines
            } else {
                Show-LineBlock -Label ($label + " (updated)") -Lines (Get-LastLines -Lines $lines -MaxLines 20)
            }
        } elseif ($count -gt 0) {
            Show-LineBlock -Label ($label + " (updated)") -Lines (Get-LastLines -Lines $lines -MaxLines 20)
        } else {
            Write-Step ($label + " was cleared or removed.") "Yellow"
        }

        $state[$path] = @{ Count = $count; Timestamp = $timestamp }
    }

    # ---- Check study directory for file changes ----
    $studyChanges = @(Check-StudyFileChanges)
    if (@($studyChanges).Count -gt 0) {
        $anyChange = $true
        Write-Host ""
        Write-Host (">>> STUDY FILE ACTIVITY (" + @($studyChanges).Count + " changes)") -ForegroundColor Yellow
        foreach ($ch in @($studyChanges)) {
            $icon = if ($ch.Type -eq "NEW") { "+" } else { "~" }
            $color = if ($ch.Type -eq "NEW") { "Green" } else { "Cyan" }
            $sizeStr = Format-FileSize $ch.Size
            Write-Host ("  " + $icon + " " + $ch.RelPath + " (" + $sizeStr + ")") -ForegroundColor $color
        }
    }

    # ---- Heartbeat (periodic status even when nothing changes) ----
    $now = Get-Date
    $elapsed = ($now - $lastHeartbeat).TotalSeconds
    if ($elapsed -ge $HeartbeatSeconds) {
        $lastHeartbeat = $now
        if (-not $anyChange) {
            $totalFiles = @($studyFileState.Keys).Count
            $ts = $now.ToString("HH:mm:ss")
            Write-Host ("[$ts] ... monitoring ($totalFiles files tracked, no changes)") -ForegroundColor DarkGray
        }
    }

    if (Test-Path -LiteralPath $PolishDoneFile) {
        Write-Host ""
        Write-Step "POLISH_COMPLETE.md detected. Study appears complete." "Green"
        break
    }

    Start-Sleep -Seconds $PollSeconds
}

Write-Host ""
Write-Step "Live monitor exiting." "Green"
