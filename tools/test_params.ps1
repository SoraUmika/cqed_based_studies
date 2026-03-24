[CmdletBinding()]
param(
    [string]$Action = 'status',
    [Parameter(Mandatory=$true)]
    [string]$StudyName,
    [string]$StudyGoal = ''
)
Write-Output "Action=$Action StudyName=$StudyName"
