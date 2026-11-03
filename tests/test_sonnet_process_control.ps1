# Dependency-free tests. Process queries and termination are mocked inside
# the imported module. This script never terminates a real process.
$ErrorActionPreference='Stop'
$module=Import-Module (Join-Path $PSScriptRoot '../tools/SonnetProcessControl.psm1') -Force -PassThru
& $module {
    $script:created=[DateTime]::Parse('2026-10-21T21:13:17.2960000Z').ToUniversalTime()
    $script:current=[pscustomobject]@{Name='em.exe';CreationDate=$script:created;
        ParentProcessId=10;SessionId=12;ProcessId=123;KernelModeTime=0;UserModeTime=10000000;WorkingSetSize=1MB}
    $script:stopCalls=0; $script:deny=$false; $script:handleTimeOffset=0; $script:nullHandle=$false
    function script:Get-CimInstance { param($ClassName,$Filter) return $script:current }
    function script:Get-Process {
        param($Id,$ErrorAction)
        if ($script:deny) { throw [ComponentModel.Win32Exception]::new(5) }
        $value=[pscustomobject]@{Handle=1;ProcessName='em';StartTime=$script:created.AddSeconds($script:handleTimeOffset)}
        if ($script:nullHandle) { $value.Handle=0; $value.StartTime=$null }
        $value | Add-Member ScriptMethod Dispose {}
        return $value
    }
    function script:Stop-Process {
        param($InputObject,[switch]$Force,$ErrorAction)
        $script:stopCalls++
    }
}
function Assert-Equal($Actual,$Expected,$Label) {
    if ($Actual -ne $Expected) { throw "$Label : expected $Expected, got $Actual" }
    Write-Output "PASS: $Label"
}
$target=[pscustomobject]@{Name='em.exe';ProcessId=123;ParentProcessId=10;SessionId=12;
    CreationUtc='2026-10-21T21:13:17.2960000Z'}
Assert-Equal (Stop-VerifiedSonnetEm $target -Preview).Status 'WouldStop' 'Preview does not terminate'
Assert-Equal (& $module {$script:stopCalls}) 0 'No stop call after preview'
& $module {$script:current.Name='python.exe'}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'IdentityChanged' 'Non-em process rejected'
& $module {$script:current.Name='em.exe';$script:current.CreationDate=$script:created.AddSeconds(1)}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'IdentityChanged' 'Reused PID rejected'
& $module {$script:current.CreationDate=$script:created;$script:handleTimeOffset=1}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'IdentityChanged' 'Handle identity rechecked'
& $module {$script:handleTimeOffset=0;$script:deny=$true}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'AccessDenied' 'Windows access denial recognized'
Assert-Equal (& $module {$script:stopCalls}) 0 'Rejected requests never terminate'
& $module {$script:deny=$false;$script:nullHandle=$true}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'AccessDenied' 'Null or zero handle treated as access denial'
& $module {$script:nullHandle=$false}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'Stopped' 'Verified em termination dispatched'
Assert-Equal (& $module {$script:stopCalls}) 1 'Exactly one mocked process termination'
& $module {$script:current=$null}
Assert-Equal (Stop-VerifiedSonnetEm $target).Status 'AlreadyExited' 'Exited process handled'
Write-Output 'All process-control checks passed. No real process was terminated.'
