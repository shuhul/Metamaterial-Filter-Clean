Set-StrictMode -Version Latest
function Invoke-SonnetEmBroker {
    [CmdletBinding()]
    param([Parameter(Mandatory)]$Request)
    $sid=[Security.Principal.WindowsIdentity]::GetCurrent().User.Value
    $pipe=[IO.Pipes.NamedPipeClientStream]::new('.',('SonnetEmControl_'+$sid),
        [IO.Pipes.PipeDirection]::InOut,[IO.Pipes.PipeOptions]::None,
        [Security.Principal.TokenImpersonationLevel]::Identification)
    try {
        $pipe.Connect(2000)
        # A normal user cannot impersonate a SYSTEM-owned pipe server.
        $owner=$pipe.GetAccessControl().GetOwner([Security.Principal.SecurityIdentifier]).Value
        if ($owner -ne 'S-1-5-18') { throw 'Refusing a broker pipe not owned by Windows SYSTEM.' }
        $json=ConvertTo-Json -InputObject $Request -Depth 8 -Compress
        $bytes=[Text.Encoding]::UTF8.GetBytes($json)
        if ($bytes.Length -gt 8192) { throw 'Request exceeds broker protocol limit.' }
        $writer=[IO.BinaryWriter]::new($pipe)
        $writer.Write([int]$bytes.Length);$writer.Write($bytes);$writer.Flush()
        $reader=[IO.BinaryReader]::new($pipe)
        $length=$reader.ReadInt32()
        if ($length -le 0 -or $length -gt 16384) { throw 'Invalid broker response length.' }
        $response=$reader.ReadBytes($length)
        if ($response.Length -ne $length) { throw 'Incomplete broker response.' }
        [Text.Encoding]::UTF8.GetString($response) | ConvertFrom-Json
    } finally { $pipe.Dispose() }
}
Export-ModuleMember -Function Invoke-SonnetEmBroker
