# /// codebase-meta
# exempt: ships inside the student kit; checks a copy of the kit on Windows with nothing installed
# ///
# Checks an unpacked course kit on Windows with Windows' own tools, for a reader without Python:
# every file against MANIFEST.json, the longest full path against Windows' limit, START.cmd's
# lookup of the page beside it, and every page START.html links, loaded by headless Edge from disk
# with every network lookup refused. Every expected value comes from the kit's own MANIFEST.json
# and START.html. Windows PowerShell 5.1; run from cmd or PowerShell in the kit's folder:
#   powershell -NoProfile -ExecutionPolicy Bypass -File check_kit.ps1
# or name another kit folder: ... -File check_kit.ps1 C:\path\to\kit
# Results also go to %TEMP%\check_kit.log, never into the kit. Exit code 0 only when every check
# met its expectation. scripts/python/check_kit.py makes the same file and link checks anywhere.

param([string]$Kit = $PSScriptRoot)

$ErrorActionPreference = 'Continue'
$Kit = (Resolve-Path -LiteralPath $Kit).Path
$log = Join-Path $env:TEMP 'check_kit.log'
Set-Content -Path $log -Value "check_kit $(Get-Date -Format s) $Kit" -Encoding ASCII
$bad = 0

function Say([string]$line) {
    Write-Output $line
    Add-Content -Path $log -Value $line -Encoding ASCII
}

$manifestFile = Join-Path $Kit 'MANIFEST.json'
if (-not (Test-Path -LiteralPath $manifestFile)) { Say "RESULT FAIL: no MANIFEST.json in $Kit"; exit 1 }
$manifest = Get-Content -LiteralPath $manifestFile -Raw -Encoding UTF8 | ConvertFrom-Json
$listed = @($manifest.entries)
Say ("RESULT kit {0}, built {1} from commit {2}; Windows {3}; PowerShell {4}" -f $Kit, $manifest.built, $manifest.commit, [Environment]::OSVersion.Version, $PSVersionTable.PSVersion)

# 1. Files: every listed file present with its size and SHA-256, and nothing unlisted.
$present = @(Get-ChildItem -LiteralPath $Kit -Recurse -File)
$names = @{}
foreach ($f in $present) { $names[$f.FullName.Substring($Kit.Length + 1).Replace('\', '/')] = $true }
$changed = 0; $missing = 0
foreach ($entry in $listed) {
    $path = Join-Path $Kit ($entry.path -replace '/', '\')
    if (-not (Test-Path -LiteralPath $path)) { $missing++; Say "  missing $($entry.path)" }
    elseif ((Get-Item -LiteralPath $path).Length -ne $entry.bytes -or
            (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLower() -ne $entry.sha256) { $changed++; Say "  changed $($entry.path)" }
    $names.Remove($entry.path)
}
$names.Remove('MANIFEST.json')
foreach ($extra in $names.Keys) { Say "  not listed $extra" }
if ($changed -gt 0 -or $missing -gt 0 -or $names.Count -gt 0) { $bad++ }
Say ("RESULT files: {0} listed; {1} changed; {2} missing; {3} not listed; expect {0}, 0, 0, 0" -f $listed.Count, $changed, $missing, $names.Count)

# 2. The longest full path where the kit sits now, against Windows' MAX_PATH of 260 with its NUL.
$longest = ($present | ForEach-Object { $_.FullName.Length } | Measure-Object -Maximum).Maximum
if ($longest -gt 259) { $bad++ }
Say ("RESULT longest full path {0} characters; expect at most 259" -f $longest)

# 3. START.cmd finds START.html beside it: a copy with "start" swapped for "echo", run from
#    another folder so the working folder cannot help it, removed even if the run fails.
$echoCmd = Join-Path $Kit 'check_kit_start_echo.cmd'
$wantStart = Join-Path $Kit 'START.html'
try {
    (Get-Content -LiteralPath (Join-Path $Kit 'START.cmd') -Raw) -replace 'start ""', 'echo' |
        Set-Content -LiteralPath $echoCmd -Encoding ASCII -NoNewline
    Push-Location $env:TEMP
    $named = (& cmd.exe /c $echoCmd | Out-String).Trim().Trim('"')
} finally {
    Pop-Location
    Remove-Item -LiteralPath $echoCmd -ErrorAction SilentlyContinue
}
if ($named -ne $wantStart -or -not (Test-Path -LiteralPath $named)) { $bad++ }
Say ("RESULT START.cmd opens {0}; exists {1}; expect {2}, True" -f $named, (Test-Path -LiteralPath $named), $wantStart)

# 4. Every page START.html links, and START.html itself, in headless Edge with every network
#    lookup refused. Each launch is bounded: --timeout caps the page, WaitForExit the process.
#    A page counts as loaded only when Edge shows its own title, since Edge's error page has one
#    too; a figure counts once, by the container Plotly adds only when its code has run.
$edge = Join-Path ${env:ProgramFiles(x86)} 'Microsoft\Edge\Application\msedge.exe'
if (-not (Test-Path -LiteralPath $edge)) {
    Say "RESULT SKIPPED pages: Edge is not at $edge; the file checks above still stand"
} else {
    $edgeData = Join-Path $env:TEMP 'check_kit_edge'
    $links = @([regex]::Matches((Get-Content -LiteralPath $wantStart -Raw -Encoding UTF8), 'href="([^"#:]+)"') | ForEach-Object { $_.Groups[1].Value })
    $noPage = 0; $short = 0; $timedOut = 0; $calls = 0; $drawnAll = 0; $i = 0
    foreach ($link in @('START.html') + $links) {
        $i++
        $file = Join-Path $Kit ($link -replace '/', '\')
        $source = Get-Content -LiteralPath $file -Raw -Encoding UTF8
        $want = [regex]::Match($source, '<title>([^<]*)</title>').Groups[1].Value
        $plotCalls = ([regex]::Matches($source, 'Plotly\.newPlot\(')).Count
        $out = Join-Path $env:TEMP ("check_kit_dom_{0}.html" -f $i)
        $arguments = "--headless=new --disable-gpu --no-first-run `"--user-data-dir=$edgeData`" `"--host-resolver-rules=MAP * ~NOTFOUND`" --virtual-time-budget=10000 --timeout=30000 --dump-dom $(([System.Uri]$file).AbsoluteUri)"
        $p = Start-Process -FilePath $edge -ArgumentList $arguments -PassThru -NoNewWindow -RedirectStandardOutput $out -RedirectStandardError "$out.err"
        $null = $p.Handle  # PowerShell 5.1 reports ExitCode only once the handle has been read
        if (-not $p.WaitForExit(120000)) { Stop-Process -Id $p.Id -Force; $timedOut++; $code = 'timeout' } else { $code = $p.ExitCode }
        $dom = ''
        if (Test-Path -LiteralPath $out) { $dom = Get-Content -LiteralPath $out -Raw -Encoding UTF8 }
        if ($dom -eq $null) { $dom = '' }
        $title = [regex]::Match($dom, '<title>([^<]*)</title>').Groups[1].Value
        if ($title -ne $want) { $noPage++ }
        $drawn = ([regex]::Matches($dom, 'class="plot-container plotly"')).Count
        if ($drawn -lt $plotCalls) { $short++ }
        $calls += $plotCalls; $drawnAll += $drawn
        Say ("  drawn {0} of {1}  exit {2}  {3} bytes  {4}  [{5}]" -f $drawn, $plotCalls, $code, $dom.Length, $title, $link)
    }
    if ($noPage -gt 0 -or $short -gt 0 -or $timedOut -gt 0) { $bad++ }
    Say ("RESULT pages {0}; without their own title {1}; figures drawn {2} of {3}; timed out {4}; expect {0}, 0, {3} of {3}, 0" -f (1 + $links.Count), $noPage, $drawnAll, $calls, $timedOut)
}
Say ("RESULT done with {0} failed checks; expect 0. The log is {1}" -f $bad, $log)
exit [int]($bad -gt 0)
