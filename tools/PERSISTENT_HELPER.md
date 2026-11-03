# Independent Sonnet EM cleanup

One administrator approval installs a narrow helper. After installation,
the assistant can stop selected Sonnet EM processes in other Windows
sessions without additional UAC prompts. No administrator password is stored.

**Installed and verified on 21 October 2026 for `DESKTOP-0ST43U4\shuhul`**
(SID ending `1064`). The SYSTEM helper rejected non-EM and stale-identity
requests, then stopped a suspended real Sonnet executable without UAC or
running a simulation. Its sampled idle CPU consumption was zero.

Evidence is in `investigations/cpu_overload_2026-10-21/persistent_helper/`.

Install once:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\setup_sonnet_em_helper.ps1
```

Check readiness, with no elevation:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1 -Action Status
```

List and stop a selected current PID, with no interactive fallback:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\sonnet_process_control.ps1 -Action Stop -ProcessIds 12345 -NoElevation
```

The ordinary tool first attempts its own process rights, then uses the helper
on access denial. `-NoElevation` permits the installed helper but prevents a
new UAC prompt if the helper is unavailable. Preview remains non-destructive.

### Agent execution identity

The Codex default shell runs as `codexsandboxoffline` (SID ending `1067`),
which is intentionally not granted broker access. For helper requests, use
`exec_command` with `sandbox_permissions="require_escalated"`; on this host
that executes as the ordinary `shuhul` account. This is host-user execution,
not Windows administrator elevation, and does not trigger UAC. Do not install
another task for the sandbox identity. Regular user PowerShell already uses
the correct account. The initial default-sandbox status timeout was resolved
by using the intended user identity; no additional privilege was granted.

## Installed components and limits

- Protected files: `C:\Program Files\SonnetEmControl\<requesting-user-SID>\`.
  SYSTEM and Administrators can modify them; the requesting user can read
  and execute them. Workspace edits cannot change installed privileged code.
- Scheduled task: `SonnetEmControl_<SID>`, running as SYSTEM, starting now
  and on system boot. Task Scheduler can restart it after failures. It waits
  on a named pipe when idle; there is no process-killing polling loop.
- Pipe access: only the requesting Windows SID, SYSTEM and Administrators.
  The client verifies that the pipe is owned by SYSTEM before sending data.
- Protocol: Status, List and Stop only. No shell commands, file paths to
  execute, arbitrary executable targets, process trees, or GUI termination.
- Stop validates the actual process image against the installed Sonnet 19.52
  `bin\em.exe` path, checks the binary hash pinned at installation, and
  checks the process creation time/session through a retained process handle.
  Renaming another executable to em.exe is insufficient. Reused PIDs fail.
- This version pins `C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe`.
  Reinstall/update the helper after a Sonnet installation change.
- Protected audit: `audit.jsonl` in the installation directory, rotated at
  approximately 5 MiB with one previous file. Workspace operation logs remain
  in `output_data/sonnet_process_control/`.

The helper is available for explicit requests; it does not independently kill
every busy job. High host CPU alone does not prove a simulation is broken.
The user's standing authorization allows the assistant to make targeted EM
cleanup requests while investigating stalled/competing runs.

## Remove or update

Remove the task and its exact protected installation directory (requires UAC):

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\tools\setup_sonnet_em_helper.ps1 -Action Remove
```

Run Install again to replace the protected code/configuration. The installer
does not change global Sonnet settings or grant the ordinary assistant process
general administrator rights. A full reinstall of Windows or policy changes
may require setup again; routine stop requests do not.
