"""Root-invoked tooling regressions: python3 -m unittest discover -s test."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SweepTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'scripts').mkdir()
        for name in ('run_cache_tier_sweep.sh', '_manifest_to_bash_arrays.py'):
            shutil.copy(ROOT / 'scripts' / name, self.root / 'scripts' / name)
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.make = self.bin / 'make'
        self.make.write_text('#!/bin/sh\nprintf "%s\\n" "DWARF_RESULT,benchmark=gemm,motif=dense-linear-algebra,target=cpu,size=8,secondary_size=8,iterations=1,repeats=1,mean_ms=1,median_ms=1,valid=${FAKE_VALID:-1}"\nexit "${FAKE_STATUS:-0}"\n')
        self.make.chmod(0o755)
        self.manifest = self.root / 'tiers.json'
        self.manifest.write_text(json.dumps({'dwarfs': {'gemm': {'l1lo': {'size': 8}}}}))

    def run_sweep(self, **overrides):
        env = dict(os.environ, PATH=f'{self.bin}:{os.environ["PATH"]}', CXX='clang++',
                   DWARF_TIERS_MANIFEST=str(self.manifest), DWARF_FILTER='gemm',
                   DWARF_TIERS='l1lo', DWARF_TARGETS='cpu')
        env.update(overrides)
        return subprocess.run(['bash', str(self.root / 'scripts/run_cache_tier_sweep.sh')],
                              cwd=self.root, env=env, capture_output=True, text=True)

    def test_empirical_benchmark_runs(self):
        result = self.run_sweep()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = next(self.root.glob('results/*/summary.csv')).read_text()
        self.assertIn('l1lo,gemm,', summary)

    def test_failed_command_is_failure_even_with_valid_result(self):
        self.assertNotEqual(self.run_sweep(FAKE_STATUS='2').returncode, 0)

    def test_invalid_result_is_failure(self):
        self.assertNotEqual(self.run_sweep(FAKE_VALID='0').returncode, 0)

    def test_missing_result_is_failure(self):
        self.make.write_text('#!/bin/sh\nexit 0\n')
        self.assertNotEqual(self.run_sweep().returncode, 0)

    def test_invalid_selections_fail(self):
        for overrides in ({'DWARF_FILTER': 'ctfft'}, {'DWARF_TIERS': 'typo'},
                          {'DWARF_TARGETS': 'typo'}, {'DWARF_FILTER': 'atax'}):
            with self.subTest(overrides=overrides):
                self.assertNotEqual(self.run_sweep(**overrides).returncode, 0)

    def test_malformed_manifest_fails(self):
        self.manifest.write_text('{broken')
        self.assertNotEqual(self.run_sweep().returncode, 0)

    def test_manifest_values_cannot_execute_shell(self):
        marker = self.root / 'executed'
        self.manifest.write_text(json.dumps({'dwarfs': {'gemm': {'l1lo': {
            'size': f'$(touch {marker})', 'secondary_size': 8}}}}))
        self.assertNotEqual(self.run_sweep().returncode, 0)
        self.assertFalse(marker.exists())


class SizingTests(unittest.TestCase):
    def test_full_suite_and_empirical_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / 'tiers.json'
            seed = Path(temp) / 'seed.json'
            seed.write_text(json.dumps({'dwarfs': {'gemm': {'l1lo': {'size': 8}}}}))
            result = subprocess.run(['python3', 'scripts/compute_cache_tier_sizes.py',
                                     '--devices', 'scripts/devices.example.json',
                                     '--seed', str(seed), '--out', str(out)],
                                    cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            dwarfs = json.loads(out.read_text())['dwarfs']
            self.assertEqual(len(dwarfs), 9)
            self.assertNotIn('ctfft', dwarfs)
            self.assertEqual(dwarfs['gemm']['l1lo']['size'], 8)

    def test_partial_streaming_seed(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / 'tiers.json'
            seed = Path(temp) / 'seed.json'
            expected = {'atax': {'l1lo': {'size': 8, 'secondary_size': 8}, 'l1hi': None}}
            seed.write_text(json.dumps({'dwarfs': expected}))
            result = subprocess.run(['python3', 'scripts/compute_cache_tier_sizes.py',
                                     '--devices', 'scripts/devices.example.json',
                                     '--seed', str(seed), '--out', str(out)],
                                    cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads(out.read_text())['dwarfs']['atax'], expected['atax'])


if __name__ == '__main__':
    unittest.main()
