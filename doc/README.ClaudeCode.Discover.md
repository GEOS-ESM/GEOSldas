# Claude Code Setup on Discover HPC

These notes capture the setup used for Claude Code on NASA Discover, where the
home directory quota is small and the installer otherwise writes too much under
`$HOME`.

## Why This Is Needed

The native Claude Code installer downloads to:

```text
$HOME/.claude/downloads
```

The installed binary also uses:

```text
$HOME/.cache/claude
$HOME/.local/share/claude
```

On Discover, `$HOME` can be around 1 GB quota. The first install attempt failed
with an `EDQUOT` quota error while copying from `~/.cache/claude/staging` into
`~/.local/share/claude/versions`.

## Check Quota

```bash
quota -s
```

For the setup captured here, `/discover/nobackup/amfox` pointed to the larger
personal nobackup fileset:

```bash
ls -ld /discover/nobackup/$USER
```

which returned:

```text
/discover/nobackup/amfox -> /gpfsm/dnb34/amfox
```

## Download Installer

Safer than piping directly to `bash`:

```bash
curl -fsSL https://claude.ai/install.sh -o install_claude.sh
```

Optionally inspect it:

```bash
less install_claude.sh
grep -i "install.*dir\\|prefix\\|CLAUDE" install_claude.sh
```

## Put Claude Storage on Nobackup

First check that these do not already contain important files:

```bash
ls -ld ~/.claude ~/.cache/claude ~/.local/share/claude
```

Create nobackup targets:

```bash
mkdir -p /discover/nobackup/amfox/.claude
mkdir -p /discover/nobackup/amfox/claude-cache
mkdir -p /discover/nobackup/amfox/claude-share
```

Remove failed partial Claude-only install directories if needed:

```bash
rm -rf ~/.cache/claude ~/.local/share/claude
```

Create parent directories and symlinks:

```bash
mkdir -p ~/.cache ~/.local/share
ln -s /discover/nobackup/amfox/.claude ~/.claude
ln -s /discover/nobackup/amfox/claude-cache ~/.cache/claude
ln -s /discover/nobackup/amfox/claude-share ~/.local/share/claude
```

Check:

```bash
ls -ld ~/.claude ~/.cache/claude ~/.local/share/claude
```

## Install

```bash
bash install_claude.sh
```

Do not use `--force` for the quota failure; the quota problem is real. Move the
Claude storage to nobackup first.

## After Install

Check size:

```bash
du -sh ~/.claude ~/.cache/claude ~/.local/share/claude
du -sh /discover/nobackup/amfox/.claude /discover/nobackup/amfox/claude-cache /discover/nobackup/amfox/claude-share
```

Check whether the launcher is available:

```bash
which claude
claude --version
```

If the installer prints shell integration or `PATH` instructions, apply them in
the shell startup file used on Discover.
