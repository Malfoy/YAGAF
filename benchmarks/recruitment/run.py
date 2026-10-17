#!/usr/bin/env python3
"""Compare recruitment/MSA settings; sequence validation uses Rust helicase tools."""
import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import time
from collections import Counter
from pathlib import Path

RECRUITMENT = {
    "adjacent31": ["--k", "31", "--flank-window", "31"],
    "window50-k31": ["--k", "31", "--flank-window", "50"],
    "window100-k31": ["--k", "31", "--flank-window", "100"],
    "dense100-k31": ["--k", "31", "--flank-window", "100", "--anchor-step", "1"],
    "window100-k25": ["--k", "25", "--flank-window", "100"],
    "dense100-k25": ["--k", "25", "--flank-window", "100", "--anchor-step", "1"],
    "window100-k21": ["--k", "21", "--flank-window", "100"],
    "dense100-k21": ["--k", "21", "--flank-window", "100", "--anchor-step", "1"],
    "best100-k31": ["--k", "31", "--flank-window", "100", "--recruitment", "best-context"],
    "best100-k25": ["--k", "25", "--flank-window", "100", "--recruitment", "best-context"],
    "best100-k21": ["--k", "21", "--flank-window", "100", "--recruitment", "best-context"],
    "best50-k31": ["--k", "31", "--flank-window", "50", "--recruitment", "best-context"],
    "window100-k41": ["--k", "41", "--flank-window", "100", "--anchor-step", "10"],
}
for strategy in ["strict","best-context"]:
    for edits in [0,1,2,4]:
        RECRUITMENT[f"{strategy}100-k31-edits{edits}"]=["--k","31","--flank-window","100","--recruitment",strategy,"--max-flank-edits",str(edits)]
RECRUITMENT["best100-k21-edits1"]=["--k","21","--flank-window","100","--recruitment","best-context","--max-flank-edits","1"]

for window,k in [(50,31),(75,31),(50,25),(50,21)]:
    RECRUITMENT[f"strict{window}-k{k}-edits1"]=["--k",str(k),"--flank-window",str(window),"--recruitment","strict","--max-flank-edits","1"]

# Preserve the permissive extension budget of the original comparison variants.
for flags in RECRUITMENT.values():
    if "--max-flank-edits" not in flags:
        flags += ["--max-flank-edits","22"]
    if "--recruitment" not in flags:
        flags += ["--recruitment","strict"]
RECRUITMENT["default"] = []

MSA = {
    f"{engine}-{cells // 1000000}M": ["--msa-backend", engine, "--max-msa-cells", str(cells), "--max-reads-per-gap", "4096"]
    for engine in ["poasta", "spoars", "spoars-scalar"]
    for cells in [20000000, 200000000, 1000000000]
}
MSA.update({f"{engine}-ont-1000M": ["--msa-backend", engine, "--preset", "ont", "--max-msa-cells", "1000000000", "--max-reads-per-gap", "4096"] for engine in ["poasta", "spoars"]})
MSA.update({f"{engine}-{order}-200M": ["--msa-backend",engine,"--msa-order",order,"--max-msa-cells","200000000","--max-reads-per-gap","4096"] for engine in ["poasta","spoars"] for order in ["longest-first","median-first"]})
MSA.update({f"{engine}-10000M": ["--msa-backend",engine,"--max-msa-cells","10000000000","--max-reads-per-gap","4096"] for engine in ["poasta","spoars","spoars-scalar"]})
MSA.update({f"{engine}-{order}-10000M": ["--msa-backend",engine,"--msa-order",order,"--max-msa-cells","10000000000","--max-reads-per-gap","4096"] for engine in ["poasta","spoars"] for order in ["longest-first","median-first"]})
MSA.update({f"{engine}-ont-10000M": ["--msa-backend",engine,"--preset","ont","--max-msa-cells","10000000000","--max-reads-per-gap","4096"] for engine in ["poasta","spoars"]})

MSA.update({f"auto-{cells // 1000000}M": ["--msa-backend","auto","--max-msa-cells",str(cells),"--max-reads-per-gap","4096"] for cells in [200000000,1000000000]})



def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scaffolds", required=True, type=Path)
    p.add_argument("--reads", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--reference", type=Path, help="Known truth for masking recovery")
    p.add_argument("--matrix", choices=["recruitment", "msa"], default="recruitment")
    p.add_argument("--only", nargs="+", help="Subset of named configurations")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--binary", type=Path, default=Path("target/release/yagaf"))
    p.add_argument("extra", nargs=argparse.REMAINDER)
    a = p.parse_args()
    configs = RECRUITMENT if a.matrix == "recruitment" else MSA
    names = a.only or list(configs)
    if any(n not in configs for n in names):
        p.error(f"Unknown configuration; choose from {list(configs)}")
    extra = a.extra[1:] if a.extra[:1] == ["--"] else a.extra
    a.output_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        dest = a.output_dir / name
        if dest.exists():
            raise SystemExit(f"Output already exists: {dest}")
        dest.mkdir()
        out = dest / "corrected.fa"
        cmd = [str(a.binary.resolve()), "-s", str(a.scaffolds), "-r", str(a.reads), "-o", str(out), "-t", str(a.threads), "--min-reads", "3", "--min-quality", "50", "--timings", str(dest / "phases.tsv")] + configs[name] + extra
        (dest / "command.json").write_text(json.dumps(cmd, indent=2) + "\n")
        before = time.monotonic()
        try:
            with (dest / "stdout.log").open("w") as stdout, (dest / "run.log").open("w") as stderr:
                process = subprocess.Popen(["/usr/bin/time", "-v", "-o", str(dest / "resources.txt")] + cmd, stdout=stdout, stderr=stderr, start_new_session=True)
                code = process.wait(timeout=a.timeout)
            status = "ok" if code == 0 else f"exit_{code}"
        except KeyboardInterrupt:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            raise
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            status = "timeout"
        summary = {"configuration": name, "status": status, "elapsed_seconds": time.monotonic() - before, "binary_sha256": hashlib.sha256(a.binary.read_bytes()).hexdigest()}
        if status == "ok":
            with Path(str(out) + ".gaps.tsv").open() as f:
                rows = list(csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"))
            summary["seed_candidate_records"] = sum(int(r["candidate_records"]) for r in rows)
            header = Path(str(out) + ".gaps.tsv").open().readline()
            summary["candidate_counting"] = next((s.split("=",1)[1] for s in header.split() if s.startswith("candidate_counting=")), "validated-intervals-only")
            if summary["candidate_counting"].endswith("paired-seeds"):
                summary["unplaced_candidates"] = sum(int(r["candidate_records"])-int(r["spanning_records"])-int(r["ambiguous_pairs"]) for r in rows)
            phases = list(csv.DictReader((dest / "phases.tsv").open(), delimiter="\t"))
            summary.update(gaps=len(rows), filled=sum(r["status"] == "filled" for r in rows), supported=sum(int(r["retained_reads"]) >= 3 for r in rows), spanning_records=sum(int(r["spanning_records"]) for r in rows), ambiguous_pairs=sum(int(r["ambiguous_pairs"]) for r in rows), reasons=dict(Counter(r["reason"] for r in rows)), phases={r["phase"]: float(r["seconds"]) for r in phases})
            for line in (dest / "resources.txt").read_text().splitlines():
                if "Maximum resident set size (kbytes):" in line:
                    summary["peak_rss_kib"] = int(line.rsplit(":", 1)[1])
            verifier = Path("target/release/examples/verify_recovery" if a.reference else "target/release/examples/verify_scaffolds")
            check = [str(verifier.resolve())]
            check += ["--reference", str(a.reference), "--masked", str(a.scaffolds), "--corrected", str(out), "--report", str(out) + ".gaps.tsv"] if a.reference else [str(a.scaffolds), str(out), str(out) + ".gaps.tsv"]
            validation = subprocess.run(check, capture_output=True, text=True, check=True)
            (dest / "verification.tsv").write_text(validation.stdout)
            (dest / "verification.log").write_text(validation.stderr)
            if a.reference:
                summary["truth"] = next(csv.DictReader(validation.stdout.splitlines(), delimiter="\t"))
        (dest / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        with (a.output_dir / "summary.jsonl").open("a") as f:
            f.write(json.dumps(summary) + "\n")
        print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
