@echo off
for %%m in (*.fastq) do for /f "tokens=1,4 delims=_" %%a in ("%%m") do ren "%%m" "%%a_%%b.fastq"
