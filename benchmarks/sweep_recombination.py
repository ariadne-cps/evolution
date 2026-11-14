#!/usr/bin/env python3
"""Run a reproducible step/group sweep of benchmark_recombination.

All summary values come from the benchmark's measured samples. The geometric
ratios describe validated outward projections, not certified set inclusion.
"""
import argparse
import csv
import json
import math
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SCENARIOS = ("urgent-crossing", "permissive-jump")
WIDTH = re.compile(r"^\s*(state\[\d+\]|time|dwell|invariant\(x-y-y\^3/3\))"
                   r"\s+independent-width=([^\s]+)\s+merged-width=([^\s]+)\s+ratio=([^\s]+)")
SAMPLE = re.compile(r"\bsample=(\d+)\s+independent-final=(\d+)"
                    r"\s+recombined-final=(\d+)\s+net-speedup=([\d.eE+-]+)")
PHASE = re.compile(r"^\s*[\u00b7|]?\s*(urgent-crossing|permissive-jump)\s+"
                   r"(contract|recombine|contract\+recombine|continue-independent|"
                   r"continue-recombined|total-independent|total-recombined)\s+\d+\s+"
                   r"([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)")

def parse(text):
    data = {name: {"at_recombination": {}, "after_continuation": {},
                   "analytic_vs_independent": {}, "analytic_vs_recombined": {},
                   "phases": {}, "samples": []} for name in SCENARIOS}
    scenario = None
    section = None
    for line in text.splitlines():
        if line.startswith("scenario="):
            scenario = line.split()[0].partition("=")[2]
            section = None
        if scenario not in data:
            continue
        if "at recombination (validated" in line:
            section = "at_recombination"
        elif "after continuation (validated" in line:
            section = "after_continuation"
        elif "analytic target flow vs independent continuation" in line:
            section = "analytic_vs_independent"
        elif "analytic target flow vs recombined continuation" in line:
            section = "analytic_vs_recombined"
        width = WIDTH.search(line)
        if width and section:
            name, before, after, ratio = width.groups()
            data[scenario][section][name] = {
                "independent_width": float(before),
                "merged_width": float(after),
                "ratio": float(ratio) if ratio != "n/a" else None}
        sample = SAMPLE.search(line)
        if sample:
            _, baseline_final, merged_final, speedup = sample.groups()
            data[scenario]["samples"].append((int(baseline_final), int(merged_final), float(speedup)))
        phase = PHASE.match(line)
        if phase:
            _, phase_name, minimum, maximum, mean, mad = phase.groups()
            data[scenario]["phases"][phase_name] = float(mean)
    return data

def run(args):
    binary = Path(args.binary).expanduser().resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"benchmark executable not found: {binary}")
    destination = Path(args.output).expanduser()
    destination.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "step", "group_size", "scenario", "status", "sample_count",
        "independent_final", "recombined_final", "speedup_mean_samples",
        "independent_mean_ns", "recombined_mean_ns", "speedup_from_means",
        "x_final_ratio", "y_final_ratio", "invariant_initial_ratio",
        "invariant_final_ratio", "analytic_x_vs_merged_ratio",
        "analytic_y_vs_merged_ratio", "analytic_invariant_vs_merged_ratio",
        "output_log", "error"]
    rows = []
    for step in args.steps:
        for size in args.groups:
            command = [str(binary), "--samples", str(args.samples),
                       "--warmups", str(args.warmups),
                       "--step", str(step), "--group-size", str(size),
                       "--accuracy", str(args.accuracy)]
            stem = f"step_{step}_group_{size}".replace(".", "p")
            logfile = destination / f"{stem}.log"
            print("RUN", " ".join(command), flush=True)
            try:
                completed = subprocess.run(command, capture_output=True, text=True,
                                           timeout=args.timeout, check=False)
                output = completed.stdout + "\n" + completed.stderr
                logfile.write_text(output, encoding="utf-8")
                parsed = parse(output)
                reason = (f"exit={completed.returncode}" if completed.returncode else "")
                for scenario in SCENARIOS:
                    part = parsed[scenario]
                    sample_rows = part["samples"]
                    phases = part["phases"]
                    clean = (completed.returncode == 0
                             and len(sample_rows) == args.samples
                             and all(key in phases for key in ("total-independent", "total-recombined"))
                             and bool(part["after_continuation"])
                             and bool(part["analytic_vs_recombined"]))
                    def ratio(section, metric):
                        return part[section].get(metric, {}).get("ratio")
                    baseline = phases.get("total-independent")
                    merged = phases.get("total-recombined")
                    row = {
                        "step": step, "group_size": size, "scenario": scenario,
                        "status": "ok" if clean else "failed",
                        "sample_count": len(sample_rows),
                        "independent_final": sample_rows[-1][0] if sample_rows else "",
                        "recombined_final": sample_rows[-1][1] if sample_rows else "",
                        "speedup_mean_samples": sum(x[2] for x in sample_rows) / len(sample_rows) if sample_rows else "",
                        "independent_mean_ns": baseline or "",
                        "recombined_mean_ns": merged or "",
                        "speedup_from_means": baseline / merged if baseline and merged else "",
                        "x_final_ratio": ratio("after_continuation", "state[0]"),
                        "y_final_ratio": ratio("after_continuation", "state[1]"),
                        "invariant_initial_ratio": ratio("at_recombination", "invariant(x-y-y^3/3)"),
                        "invariant_final_ratio": ratio("after_continuation", "invariant(x-y-y^3/3)"),
                        "analytic_x_vs_merged_ratio": ratio("analytic_vs_recombined", "state[0]"),
                        "analytic_y_vs_merged_ratio": ratio("analytic_vs_recombined", "state[1]"),
                        "analytic_invariant_vs_merged_ratio": ratio("analytic_vs_recombined", "invariant(x-y-y^3/3)"),
                        "output_log": logfile.name, "error": reason if reason else ("" if clean else "incomplete benchmark output"),
                    }
                    rows.append(row)
                    print(f"  {scenario}: {row['status']} speedup={row['speedup_from_means']}", flush=True)
            except subprocess.TimeoutExpired as exc:
                logfile.write_text(f"Timed out after {args.timeout}s\n{exc}", encoding="utf-8")
                for scenario in SCENARIOS:
                    rows.append({"step": step, "group_size": size, "scenario": scenario,
                                 "status": "timeout", "output_log": logfile.name,
                                 "error": str(exc)})
    with (destination / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    (destination / "metadata.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python": sys.version,
        "executable": str(binary), "steps": args.steps, "groups": args.groups,
        "samples": args.samples, "warmups": args.warmups,
        "accuracy": args.accuracy,
        "note": "Quality metrics are diagnostics, not inclusion certificates."
    }, indent=2) + "\n", encoding="utf-8")
    print("Wrote", destination / "results.csv", flush=True)
    return 1 if any(row["status"] != "ok" for row in rows) else 0

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="build-release/benchmarks/benchmark_recombination")
    parser.add_argument("--output", default="build-release/recombination-sweep")
    parser.add_argument("--steps", nargs="+", type=float,
                        default=[1 / 16, 1 / 64, 1 / 256])
    parser.add_argument("--groups", nargs="+", type=int, default=[2, 4, 8])
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--accuracy", type=float, default=1e-5)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    if (any(s <= 0 or not math.isfinite(s) for s in args.steps)
            or any(g < 2 for g in args.groups) or args.samples <= 0
            or args.warmups < 0 or args.timeout <= 0
            or not math.isfinite(args.accuracy) or args.accuracy <= 0):
        parser.error("steps, accuracy, samples and timeout must be positive; group sizes >= 2")
    return run(args)

if __name__ == "__main__":
    sys.exit(main())
