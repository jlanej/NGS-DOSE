#!/usr/bin/env python3
"""Turn `cargo ... --message-format=json` into GitHub annotations (and keep the readable text).

    cargo clippy --message-format=json -- -D warnings | python3 .github/scripts/cargo_annotate.py

Every compiler or clippy diagnostic becomes an `::error` / `::warning` workflow command, so it
shows up on the pull request at the line it is about. The exit status is cargo's (pipefail).
"""
import json
import sys

seen = set()
for line in sys.stdin:
    try:
        m = json.loads(line)
    except ValueError:
        sys.stdout.write(line)
        continue
    if m.get("reason") != "compiler-message":
        continue
    msg = m["message"]
    if msg.get("level") not in ("warning", "error") or not msg.get("spans"):
        if msg.get("rendered"):
            print(msg["rendered"], end="")
        continue
    span = next((s for s in msg["spans"] if s.get("is_primary")), msg["spans"][0])
    code = (msg.get("code") or {}).get("code") or msg["level"]
    key = (span["file_name"], span["line_start"], code)
    if key in seen:                      # the same finding is reported once per target (bin, test)
        continue
    seen.add(key)
    text = msg["message"].replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
    # with -D warnings every finding fails the build, so every finding is an error here
    print(f"::error file={span['file_name']},line={span['line_start']},endLine={span['line_end']},title={code}::{text}")
    print(msg.get("rendered") or msg["message"], end="")
