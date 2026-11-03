<# One-time UAC installation. Only the protected EM broker runs as SYSTEM. #>
[CmdletBinding()]
param([ValidateSet('Install','Remove','Status')][string]$Action='Install',
      [switch]$Elevated,[string]$ForUserSid)
$ErrorActionPreference='Stop'
if (-not $ForUserSid) { $ForUserSid=[Security.Principal.WindowsIdentity]::GetCurrent().User.Value }
$null=[Security.Principal.SecurityIdentifier]::new($ForUserSid)
$baseFolder=Join-Path $env:ProgramFiles 'SonnetEmControl'
$installFolder=Join-Path $baseFolder $ForUserSid
$taskName='SonnetEmControl_'+$ForUserSid
$statusFile=Join-Path $installFolder 'install_status.json'
if ($Action -eq 'Status') {
    Import-Module (Join-Path $PSScriptRoot 'SonnetEmBrokerClient.psm1') -Force
    Invoke-SonnetEmBroker -Request @{Action='Status'} | ConvertTo-Json
    exit 0
}
if (-not $Elevated) {
    $scriptLiteral=$PSCommandPath.Replace("'","''")
    $encoded=[Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes(
        "& '$scriptLiteral' -Action '$Action' -Elevated -ForUserSid '$ForUserSid'"))
    Write-Output "Requesting one-time Windows administrator approval for $Action."
    $helper=Start-Process -FilePath (Join-Path $PSHOME 'powershell.exe') -Verb RunAs -WindowStyle Hidden -PassThru `
        -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-EncodedCommand',$encoded)
    if (-not $helper.WaitForExit(30000)) { Write-Output 'Installer still running; check Status after Windows approval.'; exit 0 }
    if ($helper.ExitCode -ne 0) { throw "Installer exited with code $($helper.ExitCode). Check $statusFile" }
    if (Test-Path -LiteralPath $statusFile) { Get-Content -LiteralPath $statusFile }
    else { Write-Output "$Action completed." }
    exit 0
}
$admin=([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole(
    [Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $admin) { throw 'Installer requires Windows administrator rights.' }
# Program Files prevents ordinary users from replacing the privileged code.
if (Test-Path -LiteralPath $baseFolder) {
    if ((Get-Item -LiteralPath $baseFolder).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Install base cannot be a reparse point.' }
}
if (Test-Path -LiteralPath $installFolder) {
    if ((Get-Item -LiteralPath $installFolder).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Install folder cannot be a reparse point.' }
}
if ($Action -eq 'Remove') {
    $task=Get-ScheduledTask -TaskName $taskName -ErrorAction SilentlyContinue
    if ($task) { Stop-ScheduledTask -TaskName $taskName; Unregister-ScheduledTask -TaskName $taskName -Confirm:$false }
    if (Test-Path -LiteralPath $installFolder) {
        $resolved=(Resolve-Path -LiteralPath $installFolder).Path
        $expected=[IO.Path]::GetFullPath((Join-Path $env:ProgramFiles ('SonnetEmControl/'+$ForUserSid)))
        if ($resolved -ne $expected -or -not $resolved.StartsWith($baseFolder+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) {
            throw 'Refusing removal outside the exact helper installation directory.'
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
    exit 0
}
$allowedImage='C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe'
if (-not (Test-Path -LiteralPath $allowedImage -PathType Leaf)) { throw 'Expected Sonnet 19.52 executable missing.' }
New-Item -ItemType Directory -Path $installFolder -Force | Out-Null
$acl=[Security.AccessControl.DirectorySecurity]::new()
$acl.SetAccessRuleProtection($true,$false)
$acl.SetOwner([Security.Principal.SecurityIdentifier]::new('S-1-5-32-544'))
foreach ($entry in @(@('S-1-5-18','FullControl'),@('S-1-5-32-544','FullControl'),@($ForUserSid,'ReadAndExecute'))) {
    $acl.AddAccessRule([Security.AccessControl.FileSystemAccessRule]::new(
        [Security.Principal.SecurityIdentifier]::new($entry[0]),
        [Security.AccessControl.FileSystemRights]$entry[1],
        [Security.AccessControl.InheritanceFlags]'ContainerInherit,ObjectInherit',
        [Security.AccessControl.PropagationFlags]::None,[Security.AccessControl.AccessControlType]::Allow))
}
Set-Acl -LiteralPath $installFolder -AclObject $acl
try {
    $task=Get-ScheduledTask -TaskName $taskName -ErrorAction SilentlyContinue
    if ($task) { Stop-ScheduledTask -TaskName $taskName; Start-Sleep -Milliseconds 500 }
    $source=Join-Path $installFolder 'SonnetEmBroker.cs'
    Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'SonnetEmBroker.cs') -Destination $source -Force
    $executable=Join-Path $installFolder 'SonnetEmBroker.exe'
    $compiler=Join-Path $env:WINDIR 'Microsoft.NET/Framework64/v4.0.30319/csc.exe'
    $compilerOutput=& $compiler /nologo /optimize+ /target:exe ("/out:"+$executable) /reference:System.Web.Extensions.dll $source 2>&1
    if ($LASTEXITCODE -ne 0) { throw "Compiler failed: $compilerOutput" }
    [ordered]@{UserSid=$ForUserSid;PipeName=$taskName;AllowedImage=$allowedImage;
        AllowedHash=(Get-FileHash -LiteralPath $allowedImage -Algorithm SHA256).Hash} |
        ConvertTo-Json | Set-Content -LiteralPath (Join-Path $installFolder 'config.json') -Encoding UTF8
    $taskAction=New-ScheduledTaskAction -Execute $executable -WorkingDirectory $installFolder
    $trigger=New-ScheduledTaskTrigger -AtStartup
    $principal=New-ScheduledTaskPrincipal -UserId 'SYSTEM' -LogonType ServiceAccount -RunLevel Highest
    $settings=New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
        -StartWhenAvailable -ExecutionTimeLimit ([TimeSpan]::Zero) -RestartCount 3 `
        -RestartInterval ([TimeSpan]::FromMinutes(1)) -MultipleInstances IgnoreNew
    Register-ScheduledTask -TaskName $taskName -Action $taskAction -Trigger $trigger `
        -Principal $principal -Settings $settings -Description 'User-authorized verified Sonnet EM termination only; no shell command protocol.' -Force | Out-Null
    Start-ScheduledTask -TaskName $taskName
    [ordered]@{Status='Installed';Utc=[DateTime]::UtcNow.ToString('o');TaskName=$taskName;Folder=$installFolder;UserSid=$ForUserSid} |
        ConvertTo-Json | Set-Content -LiteralPath $statusFile -Encoding UTF8
} catch {
    [ordered]@{Status='Failed';Utc=[DateTime]::UtcNow.ToString('o');Error=$_.Exception.Message} |
        ConvertTo-Json | Set-Content -LiteralPath $statusFile -Encoding UTF8
    throw
}
