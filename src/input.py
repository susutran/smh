#!/usr/bin/env python3
import socket, os, sys

print("host:", socket.gethostname())
print("cwd:", os.getcwd())
print("python:", sys.version)
print("SLURM_NTASKS:", os.environ.get("SLURM_NTASKS"))
print("SLURM_JOB_ID:", os.environ.get("SLURM_JOB_ID"))