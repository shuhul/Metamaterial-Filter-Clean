Set-StrictMode -Version Latest

function Get-SonnetEmSnapshot {
    [CmdletBinding()]
    param([int[]]$ProcessIds)
    $items = @(Get-CimInstance Win32_Process -Filter "Name='em.exe'")
    foreach ($item in $items) {
        if ($ProcessIds -and $item.ProcessId -notin $ProcessIds) { continue }
        [pscustomobject]@{
            ProcessId = [int]$item.ProcessId
            Name = [string]$item.Name
            ParentProcessId = [int]$item.ParentProcessId
            SessionId = [int]$item.SessionId
            CreationUtc = $item.CreationDate.ToUniversalTime().ToString('o')
            CpuSeconds = ([double]$item.KernelModeTime + [double]$item.UserModeTime) / 1e7
            MemoryMiB = [Math]::Round([double]$item.WorkingSetSize / 1MB, 1)
        }
    }
}

function Stop-VerifiedSonnetEm {
    [CmdletBinding()]
    param([Parameter(Mandatory)]$Target, [switch]$Preview)
    $result = [ordered]@{ProcessId=$Target.ProcessId; Utc=[DateTime]::UtcNow.ToString('o'); Status=$null; Error=$null}
    $heldProcess = $null
    try {
        if ($Target.Name -ine 'em.exe' -or [int]$Target.ProcessId -le 0) {
            throw 'Target must identify an em.exe process.'
        }
        $current = Get-CimInstance Win32_Process -Filter "ProcessId=$([int]$Target.ProcessId)"
        if ($null -eq $current) { $result.Status='AlreadyExited'; return [pscustomobject]$result }
        if ($current.Name -ine 'em.exe' -or
            $current.CreationDate.ToUniversalTime().ToString('o') -ne $Target.CreationUtc -or
            [int]$current.ParentProcessId -ne [int]$Target.ParentProcessId -or
            [int]$current.SessionId -ne [int]$Target.SessionId) {
            $result.Status='IdentityChanged'; return [pscustomobject]$result
        }
        if ($Preview) { $result.Status='WouldStop'; return [pscustomobject]$result }
        $heldProcess = Get-Process -Id ([int]$Target.ProcessId) -ErrorAction Stop
        # Opening and retaining this handle pins the original process identity.
        # A PID reused after it exits cannot silently target a different program.
        $nativeHandle = $heldProcess.Handle
        # PowerShell can expose inaccessible Process properties as null/zero
        # instead of propagating the underlying Windows access exception.
        if ($null -eq $nativeHandle -or $nativeHandle -eq [IntPtr]::Zero -or
            $null -eq $heldProcess.StartTime) {
            $result.Status='AccessDenied'
            $result.Error='Windows did not grant a process handle/start time. Administrator rights are required.'
            return [pscustomobject]$result
        }
        $expectedTime = [DateTime]::Parse($Target.CreationUtc).ToUniversalTime()
        if ($heldProcess.ProcessName -ine 'em' -or
            [Math]::Abs(($heldProcess.StartTime.ToUniversalTime()-$expectedTime).TotalMilliseconds) -ge 1) {
            $result.Status='IdentityChanged'; return [pscustomobject]$result
        }
        Stop-Process -InputObject $heldProcess -Force -ErrorAction Stop
        $result.Status='Stopped'
    } catch {
        $result.Status='Failed'
        $exception = $_.Exception
        while ($null -ne $exception) {
            if (($exception -is [ComponentModel.Win32Exception] -and $exception.NativeErrorCode -eq 5) -or
                $exception -is [UnauthorizedAccessException]) { $result.Status='AccessDenied'; break }
            $exception = $exception.InnerException
        }
        # Stop-Process also wraps access failures in a ProcessCommandException.
        if ($_.CategoryInfo.Category -eq 'PermissionDenied') { $result.Status='AccessDenied' }
        $result.Error=$_.Exception.Message
    } finally {
        if ($null -ne $heldProcess) { $heldProcess.Dispose() }
    }
    [pscustomobject]$result
}

Export-ModuleMember -Function Get-SonnetEmSnapshot, Stop-VerifiedSonnetEm
