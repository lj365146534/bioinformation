@echo off
for %%m in (*.fastq) do for /f "tokens=1,3 delims=." %%a in ("%%m") do ren "%%m" "%%a.%%b"
