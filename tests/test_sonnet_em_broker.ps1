# Safe pre-install protocol checks: never call broker Main or stop a real EM job.
[CmdletBinding()]
param([Parameter(Mandatory)][string]$BuildPath)
$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.Web.Extensions
$assembly=[Reflection.Assembly]::LoadFile((Resolve-Path -LiteralPath $BuildPath).Path)
$type=$assembly.GetType('SonnetEmBroker')
$flags=[Reflection.BindingFlags]'NonPublic,Static'
$type.GetField('AllowedImage',$flags).SetValue($null,'C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe')
$method=$type.GetMethod('Handle',$flags)
$serializer=[Web.Script.Serialization.JavaScriptSerializer]::new()
function Check-Status($Json,$Expected,$Label) {
    $request=$serializer.DeserializeObject($Json)
    $result=$method.Invoke($null,[object[]]@($request))
    if ($result['Status'] -ne $Expected) { throw "$Label : expected $Expected, got $($result['Status'])" }
    Write-Output "PASS: $Label"
}
Check-Status '{"Action":"Execute","Command":"anything"}' 'InvalidAction' 'Arbitrary commands rejected'
Check-Status '{"Action":"Stop"}' 'InvalidTarget' 'Missing identity rejected'
Check-Status '{"Action":"Stop","Target":{"Name":"python.exe"}}' 'InvalidTarget' 'Wrong executable name rejected'
$request=@{Action='Stop';Target=@{Name='em.exe';ProcessId=$PID;SessionId=7;CreationUtc='2026-01-01T00:00:00Z'}} | ConvertTo-Json -Depth 5 -Compress
Check-Status $request 'NotAllowedExecutable' 'Live PowerShell process cannot be terminated even with a forged em name'
# Verify the client rejects an unprivileged pipe impersonating the broker.
Import-Module (Join-Path $PSScriptRoot '../tools/SonnetEmBrokerClient.psm1') -Force
$sid=[Security.Principal.WindowsIdentity]::GetCurrent().User.Value
$server=[IO.Pipes.NamedPipeServerStream]::new(('SonnetEmControl_'+$sid),[IO.Pipes.PipeDirection]::InOut,
    1,[IO.Pipes.PipeTransmissionMode]::Byte,[IO.Pipes.PipeOptions]::Asynchronous)
try {
    $pending=$server.BeginWaitForConnection($null,$null)
    $refused=$false
    try { Invoke-SonnetEmBroker -Request @{Action='Status'} | Out-Null }
    catch { if ($_.Exception.Message -like '*not owned by Windows SYSTEM*') {$refused=$true} else {throw} }
    if (-not $refused) { throw 'Client accepted a non-SYSTEM pipe' }
    if ($pending.IsCompleted) {$server.EndWaitForConnection($pending)}
    Write-Output 'PASS: Non-SYSTEM pipe impersonation rejected'
} finally { $server.Dispose() }
Write-Output 'All pre-install broker checks passed. No real process was stopped.'
