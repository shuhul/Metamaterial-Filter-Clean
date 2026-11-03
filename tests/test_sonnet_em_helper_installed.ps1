<# Installed-helper integration test. No EM simulation is run.
Creates one real em.exe process suspended before its initial thread runs;
the SYSTEM broker must reject stale identity, then terminate that idle child.
#>
$ErrorActionPreference='Stop'
Import-Module (Join-Path $PSScriptRoot '../tools/SonnetEmBrokerClient.psm1') -Force
$status=Invoke-SonnetEmBroker -Request @{Action='Status'}
if ($status.Status -ne 'Ready' -or $status.Identity -ne 'NT AUTHORITY\SYSTEM') { throw 'SYSTEM helper not ready' }
Write-Output 'PASS: Installed helper is running as SYSTEM without a UAC prompt'
$invalid=Invoke-SonnetEmBroker -Request @{Action='Execute';Command='anything'}
if ($invalid.Status -ne 'InvalidAction') { throw 'Arbitrary command was not rejected' }
$self=Get-CimInstance Win32_Process -Filter "ProcessId=$PID"
$fake=@{Name='em.exe';ProcessId=$PID;SessionId=[int]$self.SessionId;CreationUtc=$self.CreationDate.ToUniversalTime().ToString('o')}
$rejected=Invoke-SonnetEmBroker -Request @{Action='Stop';Target=$fake}
if ($rejected.Status -ne 'NotAllowedExecutable') { throw 'Non-EM target was not rejected' }
Write-Output 'PASS: Privileged helper rejects arbitrary commands and a live non-EM process'
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class SuspendedSonnetTest {
 [StructLayout(LayoutKind.Sequential,CharSet=CharSet.Unicode)]
 public struct Startup {
  public int cb; public string reserved, desktop, title;
  public int x,y,xSize,ySize,xChars,yChars,fill,flags;
  public short show,reservedBytes; public IntPtr reservedPointer,input,output,error;
 }
 [StructLayout(LayoutKind.Sequential)] public struct Info {
  public IntPtr process,thread; public int processId,threadId;
 }
 [DllImport("kernel32.dll",CharSet=CharSet.Unicode,SetLastError=true)]
 public static extern bool CreateProcess(string application,string command,IntPtr pa,IntPtr ta,
   bool inherit,uint flags,IntPtr environment,string directory,ref Startup startup,out Info info);
 [DllImport("kernel32.dll")] public static extern bool TerminateProcess(IntPtr process,uint code);
 [DllImport("kernel32.dll")] public static extern uint WaitForSingleObject(IntPtr process,uint milliseconds);
 [DllImport("kernel32.dll")] public static extern bool CloseHandle(IntPtr handle);
}
'@
$start=[SuspendedSonnetTest+Startup]::new()
$start.cb=[Runtime.InteropServices.Marshal]::SizeOf($start)
$info=[SuspendedSonnetTest+Info]::new()
try {
    # CREATE_SUSPENDED | CREATE_NO_WINDOW; no Sonnet engine instructions execute.
    $created=[SuspendedSonnetTest]::CreateProcess($status.AllowedImage,('"'+$status.AllowedImage+'"'),
        [IntPtr]::Zero,[IntPtr]::Zero,$false,0x08000004,[IntPtr]::Zero,
        (Split-Path $status.AllowedImage -Parent),[ref]$start,[ref]$info)
    if (-not $created) { throw [ComponentModel.Win32Exception]::new([Runtime.InteropServices.Marshal]::GetLastWin32Error()) }
    $list=Invoke-SonnetEmBroker -Request @{Action='List'}
    $target=@($list.Processes | Where-Object ProcessId -eq $info.processId)
    if ($target.Count -ne 1) { throw 'Suspended child not found through helper' }
    $wrong=@{Name='em.exe';ProcessId=$info.processId;SessionId=$target[0].SessionId;
        CreationUtc=([DateTime]::Parse($target[0].CreationUtc).AddSeconds(1).ToString('o'))}
    $mismatch=Invoke-SonnetEmBroker -Request @{Action='Stop';Target=$wrong}
    if ($mismatch.Status -ne 'IdentityChanged' -or [SuspendedSonnetTest]::WaitForSingleObject($info.process,0) -ne 258) {
        throw 'Stale identity was not rejected without termination'
    }
    Write-Output 'PASS: Stale creation time rejected; suspended EM child remains alive'
    $result=Invoke-SonnetEmBroker -Request @{Action='Stop';Target=$target[0]}
    if ($result.Status -ne 'Stopped' -or [SuspendedSonnetTest]::WaitForSingleObject($info.process,1000) -ne 0) {
        throw 'Helper did not stop the verified suspended EM child'
    }
    Write-Output "PASS: SYSTEM helper stopped verified idle EM child $($info.processId), without UAC or an EM simulation"
} finally {
    if ($info.process -ne [IntPtr]::Zero) {
        if ([SuspendedSonnetTest]::WaitForSingleObject($info.process,0) -eq 258) {
            [SuspendedSonnetTest]::TerminateProcess($info.process,1) | Out-Null
        }
        [SuspendedSonnetTest]::CloseHandle($info.process) | Out-Null
    }
    if ($info.thread -ne [IntPtr]::Zero) { [SuspendedSonnetTest]::CloseHandle($info.thread) | Out-Null }
}
