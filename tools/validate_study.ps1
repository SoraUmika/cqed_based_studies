# validate_study.ps1
# Automated pre-review structural validator for cQED studies.
# Checks mandatory files, README sections, IMPROVEMENTS.md tags,
# artifacts, figures, and report compilation.
#
# Usage:
#   .\tools\validate_study.ps1 -StudyName <name>
#   .\tools\validate_study.ps1 -StudyPath studies\<name>
#
# Returns exit code 0 if all checks pass, 1 otherwise.

param(
    [string]$StudyName,
    [string]$StudyPath
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Resolve study path
if ($StudyPath) {
    $StudyDir = Resolve-Path $StudyPath -ErrorAction SilentlyContinue
} elseif ($StudyName) {
    $StudyDir = Join-Path $PSScriptRoot ".."
    $StudyDir = Join-Path $StudyDir "studies"
    $StudyDir = Join-Path $StudyDir $StudyName
} else {
    Write-Error "Provide -StudyName or -StudyPath"
    exit 1
}

if (-not (Test-Path $StudyDir)) {
    Write-Error "Study directory not found: $StudyDir"
    exit 1
}

$pass = 0
$fail = 0
$warn = 0

function Check-Pass {
    param([string]$Name)
    $script:pass++
    Write-Output "  [PASS] $Name"
}

function Check-Fail {
    param([string]$Name, [string]$Detail)
    $script:fail++
    Write-Output "  [FAIL] $Name -- $Detail"
}

function Check-Warn {
    param([string]$Name, [string]$Detail)
    $script:warn++
    Write-Output "  [WARN] $Name -- $Detail"
}

Write-Output ""
Write-Output ("=" * 64)
Write-Output "  Study Validation: $(Split-Path $StudyDir -Leaf)"
Write-Output ("=" * 64)
Write-Output ""

# ──────────────────────────────────────────────────────────
# 1. Mandatory Files
# ──────────────────────────────────────────────────────────
Write-Output "--- 1. Mandatory Files ---"

$mandatoryFiles = @(
    "README.md",
    "IMPROVEMENTS.md"
)

foreach ($f in $mandatoryFiles) {
    $fp = Join-Path $StudyDir $f
    if (Test-Path $fp) {
        Check-Pass "$f exists"
    } else {
        Check-Fail "$f exists" "File not found"
    }
}

# Directories
$mandatoryDirs = @("scripts", "data", "figures", "artifacts", "report")
foreach ($d in $mandatoryDirs) {
    $dp = Join-Path $StudyDir $d
    if (Test-Path $dp) {
        Check-Pass "$d/ directory exists"
    } else {
        Check-Fail "$d/ directory exists" "Directory not found"
    }
}

# Report files
$reportTex = Join-Path (Join-Path $StudyDir "report") "report.tex"
if (Test-Path $reportTex) {
    Check-Pass "report/report.tex exists"
} else {
    Check-Fail "report/report.tex exists" "File not found"
}

$reportPdf = Join-Path (Join-Path $StudyDir "report") "report.pdf"
if (Test-Path $reportPdf) {
    Check-Pass "report/report.pdf exists"
} else {
    Check-Warn "report/report.pdf exists" "PDF not compiled yet"
}

# Reproducibility notebook
$notebookPath = Join-Path (Join-Path $StudyDir "scripts") "reproducibility_notebook.ipynb"
if (Test-Path $notebookPath) {
    Check-Pass "scripts/reproducibility_notebook.ipynb exists"
} else {
    Check-Fail "scripts/reproducibility_notebook.ipynb exists" "Required by AGENTS.md Step 6"
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# 2. README.md Sections
# ──────────────────────────────────────────────────────────
Write-Output "--- 2. README.md Sections ---"

$readmePath = Join-Path $StudyDir "README.md"
if (Test-Path $readmePath) {
    $readmeContent = Get-Content $readmePath -Raw -Encoding UTF8

    $requiredSections = @(
        "## Problem Class",
        "## Motivation",
        "## Goals",
        "## Methods",
        "## Analytic Preliminary",
        "## cqed_sim Gap Analysis",
        "## Assumptions",
        "## Compute & Resource Strategy",
        "## Expected Outcomes",
        "## Known Limitations",
        "## Validation",
        "## Status"
    )

    foreach ($section in $requiredSections) {
        if ($readmeContent -match [regex]::Escape($section)) {
            Check-Pass "README has '$section'"
        } else {
            Check-Fail "README has '$section'" "Section missing"
        }
    }

    # Check validation checkboxes
    if ($readmeContent -match "\[[ x]\] Sanity checks") {
        Check-Pass "README Validation has sanity checks item"
    } else {
        Check-Warn "README Validation has sanity checks item" "Expected checkbox"
    }

    if ($readmeContent -match "\[[ x]\] Convergence") {
        Check-Pass "README Validation has convergence item"
    } else {
        Check-Warn "README Validation has convergence item" "Expected checkbox"
    }
} else {
    Check-Fail "README.md sections" "File not found"
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# 3. IMPROVEMENTS.md Structure
# ──────────────────────────────────────────────────────────
Write-Output "--- 3. IMPROVEMENTS.md ---"

$improvPath = Join-Path $StudyDir "IMPROVEMENTS.md"
if (Test-Path $improvPath) {
    $improvContent = Get-Content $improvPath -Raw -Encoding UTF8

    $improvSections = @(
        "## Critical Gaps (P1)",
        "## Recommended Improvements (P2)",
        "## Nice-to-Haves (P3)",
        "## What Was Tried and Did Not Work",
        "## Compute & Resource Notes"
    )

    foreach ($section in $improvSections) {
        if ($improvContent -match [regex]::Escape($section)) {
            Check-Pass "IMPROVEMENTS has '$section'"
        } else {
            Check-Fail "IMPROVEMENTS has '$section'" "Section missing"
        }
    }

    # Check for priority/difficulty tags
    $tagPattern = "\*\*\[P[123]\s*[\|/]\s*(LOW|MEDIUM|HIGH)\]\*\*|\[P[123]\].*\[(LOW|MEDIUM|HIGH)\]|\(P[123],\s*(LOW|MEDIUM|HIGH)\)"
    $hasItems = ($improvContent -match "^- " -or $improvContent -match "^\* ")
    if ($improvContent.Length -gt 200 -and -not ($improvContent -match $tagPattern)) {
        Check-Warn "IMPROVEMENTS has priority/difficulty tags" "Items may lack P1/P2/P3 + LOW/MEDIUM/HIGH tags"
    } else {
        Check-Pass "IMPROVEMENTS format looks reasonable"
    }
} else {
    Check-Fail "IMPROVEMENTS.md" "File not found"
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# 4. Artifacts
# ──────────────────────────────────────────────────────────
Write-Output "--- 4. Artifacts ---"

$artifactsDir = Join-Path $StudyDir "artifacts"
if (Test-Path $artifactsDir) {
    $artifacts = @(Get-ChildItem $artifactsDir -File -ErrorAction SilentlyContinue)
    if ($artifacts.Count -gt 0) {
        Check-Pass "artifacts/ has $($artifacts.Count) file(s)"

        # Check JSON artifacts for metadata
        $jsonFiles = @($artifacts | Where-Object { $_.Extension -eq ".json" })
        if ($jsonFiles.Count -gt 0) {
            foreach ($jf in $jsonFiles | Select-Object -First 3) {
                try {
                    $content = Get-Content $jf.FullName -Raw -Encoding UTF8
                    $null = $content | ConvertFrom-Json
                    Check-Pass "$($jf.Name) is valid JSON"
                } catch {
                    Check-Fail "$($jf.Name) is valid JSON" "Parse error"
                }
            }
        }
    } else {
        Check-Fail "artifacts/ has files" "Directory is empty"
    }
} else {
    Check-Fail "artifacts/ directory" "Not found"
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# 5. Figures (dual format)
# ──────────────────────────────────────────────────────────
Write-Output "--- 5. Figures ---"

$figuresDir = Join-Path $StudyDir "figures"
if (Test-Path $figuresDir) {
    $figFiles = @(Get-ChildItem $figuresDir -File -ErrorAction SilentlyContinue)
    if ($figFiles.Count -gt 0) {
        Check-Pass "figures/ has $($figFiles.Count) file(s)"

        # Check for dual format
        $pngFiles = @($figFiles | Where-Object { $_.Extension -eq ".png" })
        $pdfFiles = @($figFiles | Where-Object { $_.Extension -eq ".pdf" })

        if ($pngFiles.Count -gt 0 -and $pdfFiles.Count -gt 0) {
            Check-Pass "figures/ has both PNG and PDF formats"
        } elseif ($pngFiles.Count -gt 0 -and $pdfFiles.Count -eq 0) {
            Check-Warn "figures/ dual format" "PNG found but no PDF (vector); save both per AGENTS.md"
        } elseif ($pdfFiles.Count -gt 0 -and $pngFiles.Count -eq 0) {
            Check-Warn "figures/ dual format" "PDF found but no PNG (raster); save both per AGENTS.md"
        }

        # Check for descriptive names
        $numberedNames = @($figFiles | Where-Object { $_.BaseName -match "^fig\d+$" })
        if ($numberedNames.Count -gt 0) {
            Check-Warn "Figure naming" "Found numbered names (fig1, fig2...); use descriptive names"
        }
    } else {
        Check-Fail "figures/ has files" "Directory is empty"
    }
} else {
    Check-Fail "figures/ directory" "Not found"
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# 6. Report Checks
# ──────────────────────────────────────────────────────────
Write-Output "--- 6. Report ---"

if (Test-Path $reportTex) {
    $texContent = Get-Content $reportTex -Raw -Encoding UTF8

    # Check for mandatory sections
    $mandatorySections = @(
        "\begin{abstract}",
        "\section{Introduction}",
        "\section{Validation}",
        "\section{Discussion}",
        "\section{Conclusion}",
        "\section{Limitations"
    )

    foreach ($s in $mandatorySections) {
        $escaped = [regex]::Escape($s)
        if ($texContent -match $escaped) {
            Check-Pass "report.tex has '$s'"
        } else {
            Check-Fail "report.tex has '$s'" "Section missing"
        }
    }

    # Check for appendix
    if ($texContent -match "\\appendix") {
        Check-Pass "report.tex has appendix"
    } else {
        Check-Fail "report.tex has appendix" "Appendix is mandatory (AGENTS.md)"
    }

    # Check for filenames in main text (before appendix)
    $beforeAppendix = $texContent
    $appIdx = $texContent.IndexOf("\appendix")
    if ($appIdx -gt 0) {
        $beforeAppendix = $texContent.Substring(0, $appIdx)
    }

    $filenamePatterns = @("\.py\b", "\.json\b", "\.npz\b", "\.csv\b", "data/", "artifacts/", "scripts/")
    foreach ($pat in $filenamePatterns) {
        if ($beforeAppendix -match $pat) {
            Check-Warn "No filenames in main text" "Pattern '$pat' found before appendix"
        }
    }

    # Check for snake_case in prose (rough heuristic)
    $snakePattern = "[a-z]+_[a-z]+_[a-z]+"
    $matchCount = @([regex]::Matches($beforeAppendix, $snakePattern)).Count
    if ($matchCount -gt 5) {
        Check-Warn "No snake_case in main text" "Found ~$matchCount potential snake_case identifiers before appendix"
    }
}

Write-Output ""

# ──────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────
Write-Output ("=" * 64)
Write-Output "  RESULTS: $pass PASS, $fail FAIL, $warn WARN"
Write-Output ("=" * 64)

if ($fail -eq 0) {
    Write-Output ""
    Write-Output "  All structural checks passed. Study is ready for review."
    exit 0
} else {
    Write-Output ""
    Write-Output "  $fail check(s) failed. Fix before writing REVIEW_REQUEST.md."
    exit 1
}
