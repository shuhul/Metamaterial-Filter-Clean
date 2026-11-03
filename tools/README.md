# Sonnet process control

For independent cross-session cleanup after a one-time administrator setup,
use the [persistent helper](PERSISTENT_HELPER.md). The existing stop tool
automatically uses it on access denial before considering another UAC prompt.

Run from the project root in PowerShell. No additional packages are needed.

List EM processes:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1
```

Preview a stop, replacing the example PID with the current one:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1 -Action Stop -ProcessIds 46708 -Preview
```

Stop it:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1 -Action Stop -ProcessIds 46708
```

The stop command tries current credentials first. Access denial triggers a
**Windows UAC administrator prompt**. Approve it or provide authorized admin
credentials. An already elevated terminal normally needs no second prompt.
`-NoElevation` records access denial without requesting UAC.

Permission given to the assistant does not confer Windows administrator
rights. This tool cannot bypass UAC. Without admin approval or the owning
user's credentials, another user's process may remain inaccessible.

Only explicitly selected `em.exe` processes are eligible. It rechecks the
name, creation time, session and parent PID, then retains a process handle
and rechecks its name/start time before stopping it. Reused PIDs are rejected.
There is no kill-all default, process-tree termination or GUI termination.
A separate launcher could start another EM job after its child is stopped.

Each operation saves a snapshot and results under
`output_data/sonnet_process_control/<timestamp>_<id>/`. Elevated requests add
`request.json`, then `elevated_result.json` or `elevation_error.json`.
The per-request UAC helper exits afterward. The separately installed
[persistent helper](PERSISTENT_HELPER.md) remains available for future requests.

Test without terminating real processes:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tests\test_sonnet_process_control.ps1
```

See the [CPU investigation](../investigations/cpu_overload_2026-10-21/README.md).

The elevated stop path was verified on 21 October 2026: PID 46708 was stopped
and measured host CPU dropped from 96.32% to 3.01%. That PID is an example,
not a permanent target; always list current processes before a new request.
