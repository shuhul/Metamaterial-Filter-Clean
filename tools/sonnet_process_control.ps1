<#
List:    .\tools\sonnet_process_control.ps1
Preview: .\tools\sonnet_process_control.ps1 -Action Stop -ProcessIds 46708 -Preview
Stop:    .\tools\sonnet_process_control.ps1 -Action Stop -ProcessIds 46708

Stop requests Windows UAC elevation only if access is denied. Approving the
Windows prompt (or administrator credentials) is required. There is no bypass.
Only explicit em.exe snapshots can be stopped; no GUI, Python, or process tree.
#>
[CmdletBinding()]
param(
    [ValidateSet('List','Stop','Status')][string]$Action='List',
    [int[]]$ProcessIds,
    [switch]$Preview,
    [switch]$NoElevation,
    [string]$RequestFile
)
$ErrorActionPreference='Stop'
Import-Module (Join-Path $PSScriptRoot 'SonnetProcessControl.psm1') -Force
$workspaceRoot = Split-Path $PSScriptRoot -Parent
$auditRoot = Join-Path $workspaceRoot 'output_data/sonnet_process_control'
New-Item -ItemType Directory -Path $auditRoot -Force | Out-Null

if ($Action -eq 'Status') {
    Import-Module (Join-Path $PSScriptRoot 'SonnetEmBrokerClient.psm1') -Force
    Invoke-SonnetEmBroker -Request @{Action='Status'} | ConvertTo-Json
    exit 0
}

function Save-Json($Path, $Value) {
    ConvertTo-Json -InputObject $Value -Depth 8 | Set-Content -LiteralPath $Path -Encoding UTF8
}

if ($RequestFile) {
    # Elevated worker accepts only a saved, local request, never shell commands.
    $resolvedRequest = (Resolve-Path -LiteralPath $RequestFile).Path
    $resolvedAuditRoot = (Resolve-Path -LiteralPath $auditRoot).Path + [IO.Path]::DirectorySeparatorChar
    if (-not $resolvedRequest.StartsWith($resolvedAuditRoot,[StringComparison]::OrdinalIgnoreCase)) {
        throw 'Request must be inside this workspace process-control audit folder.'
    }
    $isAdmin = ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole(
        [Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $isAdmin) { throw 'Elevated worker requires Windows administrator rights.' }
    $request = Get-Content -LiteralPath $resolvedRequest -Raw | ConvertFrom-Json
    $results = @($request.Targets | ForEach-Object { Stop-VerifiedSonnetEm -Target $_ })
    $resultFile = Join-Path (Split-Path $resolvedRequest -Parent) 'elevated_result.json'
    Save-Json $resultFile $results
    $results | Format-Table ProcessId,Status,Error -AutoSize
    exit 0
}

$runFolder = Join-Path $auditRoot (([DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffZ'))+'_'+[Guid]::NewGuid().ToString('N').Substring(0,8))
New-Item -ItemType Directory -Path $runFolder | Out-Null
$targets = @(Get-SonnetEmSnapshot -ProcessIds $ProcessIds)
Save-Json (Join-Path $runFolder 'snapshot.json') $targets
if ($Action -eq 'List') {
    $targets | Format-Table ProcessId,ParentProcessId,SessionId,MemoryMiB,CpuSeconds,CreationUtc -AutoSize
    Write-Output "Audit folder: $runFolder"
    exit 0
}
if (-not $ProcessIds -or @($ProcessIds | Where-Object { $_ -le 0 }).Count) {
    throw 'Stop requires explicit positive -ProcessIds. Use List first.'
}
$results = @($targets | ForEach-Object { Stop-VerifiedSonnetEm -Target $_ -Preview:$Preview })
foreach ($requestedId in $ProcessIds) {
    if ($requestedId -notin @($targets | ForEach-Object { $_.ProcessId })) {
        $results += [pscustomobject]@{ProcessId=$requestedId;Utc=[DateTime]::UtcNow.ToString('o');Status='NotAnEmProcessOrExited';Error=$null}
    }
}
Save-Json (Join-Path $runFolder 'result.json') $results
$results | Format-Table ProcessId,Status,Error -AutoSize
$denied = @($results | Where-Object Status -eq 'AccessDenied')
if ($denied.Count -and -not $Preview) {
    # Prefer the one-time installed narrow SYSTEM helper. No UAC per request.
    Import-Module (Join-Path $PSScriptRoot 'SonnetEmBrokerClient.psm1') -Force
    foreach ($deniedItem in $denied) {
        try {
            $target=@($targets | Where-Object ProcessId -eq $deniedItem.ProcessId)[0]
            $brokerResult=Invoke-SonnetEmBroker -Request @{Action='Stop';Target=$target}
            $deniedItem.Status=$brokerResult.Status
            $deniedItem.Error=if ($brokerResult.PSObject.Properties['Error']) { $brokerResult.Error } else { $null }
        } catch {
            # Missing/unavailable helper retains the existing UAC fallback.
            $deniedItem.Error=$deniedItem.Error+' Broker unavailable: '+$_.Exception.Message
        }
    }
    Save-Json (Join-Path $runFolder 'result.json') $results
    $results | Format-Table ProcessId,Status,Error -AutoSize
    $denied = @($results | Where-Object Status -eq 'AccessDenied')
}
if ($denied.Count -and -not $Preview -and -not $NoElevation) {
    $requestPath = Join-Path $runFolder 'request.json'
    $deniedIds = @($denied | ForEach-Object { $_.ProcessId })
    Save-Json $requestPath ([ordered]@{CreatedUtc=[DateTime]::UtcNow.ToString('o');
        Targets=@($targets | Where-Object { $_.ProcessId -in $deniedIds })})
    # Single-quoted literal paths + EncodedCommand preserve spaces/apostrophes.
    $literalScript = $PSCommandPath.Replace("'","''")
    $literalRequest = $requestPath.Replace("'","''")
    $encoded = [Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes(
        "& '$literalScript' -RequestFile '$literalRequest'"))
    Write-Output 'Windows administrator approval is required; requesting the UAC prompt.'
    try {
        $helper = Start-Process -FilePath (Join-Path $PSHOME 'powershell.exe') -Verb RunAs -WindowStyle Hidden -PassThru `
            -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-EncodedCommand',$encoded)
        if ($helper.WaitForExit(30000)) {
            $elevatedResult = Join-Path $runFolder 'elevated_result.json'
            if (Test-Path -LiteralPath $elevatedResult) {
                Get-Content -LiteralPath $elevatedResult -Raw | ConvertFrom-Json | Format-Table ProcessId,Status,Error -AutoSize
            } else { Write-Warning 'Helper exited without a result; inspect this audit folder.' }
        } else { Write-Output 'Elevated helper still running; inspect elevated_result.json when it finishes.' }
    } catch {
        Save-Json (Join-Path $runFolder 'elevation_error.json') ([ordered]@{
            Utc=[DateTime]::UtcNow.ToString('o');Error=$_.Exception.Message})
        Write-Warning "Windows elevation did not complete: $($_.Exception.Message)"
    }
}
Write-Output "Audit folder: $runFolder"
