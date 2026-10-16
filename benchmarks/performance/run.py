#!/usr/bin/env python3
"""Run serial, alternating baseline/current comparisons on prepared recovery data."""
import argparse
import csv
import hashlib
import pathlib
import statistics
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('--baseline', type=pathlib.Path, required=True)
p.add_argument('--binary', type=pathlib.Path, required=True)
p.add_argument('--data', type=pathlib.Path, required=True)
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repeats', type=int, default=7)
p.add_argument('--threads', type=int, default=32)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
rows = []
phase_rows = []

def digest(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()

def report(path):
    return [line for line in path.read_text().splitlines() if not line.startswith('#')]

for strain in ('mg1655', 'w3110'):
    for masks in (1000, 10000):
        dataset = a.data / f'{strain}-50x-mask{masks}-seed42'
        expected = None
        for repetition in range(a.repeats):
            order = ('baseline', 'optimized') if repetition % 2 == 0 else ('optimized', 'baseline')
            for version in order:
                run = a.output / f'{strain}-{masks}-{version}-{repetition + 1}'
                run.mkdir()
                binary = a.baseline if version == 'baseline' else a.binary
                command = [str(binary.resolve()), '--scaffolds', str((dataset / 'masked.fa').resolve()),
                    '--reads', str((dataset / 'reads.fa').resolve()), '--output', str(run / 'corrected.fa'),
                    '--report', str(run / 'gaps.tsv'), '--timings', str(run / 'phases.tsv'),
                    '--preset', 'hifi', '--k', '31', '--min-reads', '5', '--min-quality', '100',
                    '--threads', str(a.threads), '--msa-threads', str(a.threads)]
                start = time.perf_counter()
                with (run / 'stderr.txt').open('w') as stderr:
                    subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(run / 'rss.txt'), *command],
                        stdout=subprocess.DEVNULL, stderr=stderr, check=True)
                elapsed = time.perf_counter() - start
                signature = (digest(run / 'corrected.fa'), report(run / 'gaps.tsv'))
                if expected is None:
                    expected = signature
                assert signature == expected, f'Results changed: {run}'
                row = dict(strain=strain, masked_bases=masks, version=version, repetition=repetition + 1,
                    threads=a.threads, msa_threads=a.threads, process_seconds=f'{elapsed:.9f}',
                    peak_rss_kib=(run / 'rss.txt').read_text().strip(), fasta_sha256=signature[0])
                rows.append(row)
                for phase in csv.DictReader((run / 'phases.tsv').open(), delimiter='\t'):
                    phase_rows.append(dict(strain=strain, masked_bases=masks, version=version,
                        repetition=repetition + 1, **phase))
                print(strain, masks, version, repetition + 1, f'{elapsed:.4f}s', flush=True)
                # Retain small reports/timings; do not duplicate large corrected scaffolds.
                (run / 'corrected.fa').unlink()

for filename, records in [('runs.tsv', rows), ('phases.tsv', phase_rows)]:
    with (a.output / filename).open('w') as output:
        writer = csv.DictWriter(output, fieldnames=records[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
for strain in ('mg1655', 'w3110'):
    for masks in (1000, 10000):
        for version in ('baseline', 'optimized'):
            values = [float(r['process_seconds']) for r in rows
                if r['strain'] == strain and r['masked_bases'] == masks and r['version'] == version]
            print(strain, masks, version, 'median', statistics.median(values))
