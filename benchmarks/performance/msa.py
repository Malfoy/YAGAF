#!/usr/bin/env python3
"""Compare native MSA engines on the same deterministic reference-derived spans."""
import argparse
import csv
import pathlib
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--binary', type=pathlib.Path, required=True, help='release examples/msa_bench')
p.add_argument('--reference', type=pathlib.Path, required=True)
p.add_argument('--output', type=pathlib.Path, required=True)
p.add_argument('--repeats', type=int, default=5)
a = p.parse_args()
a.output.mkdir(parents=True, exist_ok=False)
rows = []
for preset, error in [('hifi', 0.001), ('ont', 0.05)]:
    for length, depth in [(128, 50), (500, 50), (1000, 16)]:
        scalar_metrics = None
        for backend in ['poasta', 'spoars', 'spoars-scalar']:
            stem = a.output / f'{preset}-{length}-{backend}'
            command = [str(a.binary.resolve()), '--reference', str(a.reference.resolve()),
                '--backend', backend, '--preset', preset, '--error-rate', str(error),
                '--output', str(stem.with_suffix('.tsv')), '--length', str(length),
                '--depth', str(depth), '--families', '32', '--threads', '32',
                '--repeats', str(a.repeats)]
            with stem.with_suffix('.stderr').open('w') as stderr:
                subprocess.run(['/usr/bin/time', '-f', '%M', '-o', str(stem.with_suffix('.rss')), *command],
                    stdout=subprocess.DEVNULL, stderr=stderr, check=True)
            result = list(csv.DictReader(stem.with_suffix('.tsv').open(), delimiter='\t'))
            metrics = [[r[k] for k in ('exact_consensuses', 'rejected_consensuses', 'unanimous_columns', 'msa_columns')]
                for r in result]
            assert all(m == metrics[0] for m in metrics), 'Nondeterministic MSA metrics'
            if backend == 'spoars':
                scalar_metrics = metrics
            if backend == 'spoars-scalar':
                assert metrics == scalar_metrics, 'SIMD/scalar results differ'
            for row in result:
                rows.append(dict(**row, peak_rss_kib=stem.with_suffix('.rss').read_text().strip()))
            print(preset, length, backend, [r['seconds'] for r in result], metrics[0], flush=True)
with (a.output / 'runs.tsv').open('w') as output:
    writer = csv.DictWriter(output, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
