#!/usr/bin/env python3
"""CLI and input regressions for already-built benchmark binaries."""
import argparse
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--device', default='cpu')
    parser.add_argument('--build-dir', default='build')
    args = parser.parse_args()
    paths = (ROOT / 'common/benchmarks.mk').read_text().replace(
        'BENCHMARK_DIRS :=', '').replace('\\', '').split()
    binaries = {Path(path).name: ROOT / path / args.build_dir /
                f'{Path(path).name}_{args.device}' for path in paths}
    checks = 0

    def run(name, flags, success=True):
        nonlocal checks
        result = subprocess.run([str(binaries[name]), *flags], cwd=ROOT,
                                capture_output=True, text=True, timeout=60)
        assert (result.returncode == 0) == success, (
            name, flags, result.returncode, result.stdout, result.stderr)
        checks += 1
        return result.stdout

    for name in binaries:
        assert 'Options:' in run(name, ['--help'])
        for flags in (['--typo'], ['--size'], ['--size', 'abc'],
                      ['--size', '-1'], ['--repeats', '0'],
                      ['--iterations', '1junk']):
            run(name, flags, success=False)

    small = ['--size', '8', '--secondary-size', '8', '--warmup', '0', '--repeats', '1']
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        empty = tmp / 'empty.txt'
        empty.write_text('')
        events, ivl, episodes = [tmp / name for name in ('events.csv', 'ivl.txt', 'episodes.txt')]
        events.write_text('one,two,three\nfour\n')
        ivl.write_text('a b\n')
        episodes.write_text('c d\n')
        flags = small + ['--input', str(events), '--input2', str(ivl), '--input3', str(episodes)]
        output = run('tdm', flags)
        result = next(line for line in output.splitlines() if line.startswith('DWARF_RESULT,'))
        fields = dict(item.split('=', 1) for item in result.split(',')[1:])
        assert fields['size'] == '4' and fields['secondary_size'] == '8', result
        run('tdm', small + ['--input', str(empty), '--input2', str(empty),
                           '--input3', str(empty)], success=False)
        run('tdm', small + ['--input', str(events), '--input2', str(tmp / 'missing'),
                           '--input3', str(episodes)], success=False)
        # Source 2 is isolated; a valid traversal visits only that node.
        graph = tmp / 'graph.txt'
        graph.write_text('3\n0 1\n1 0\n1 0\n2\n1\n1 1\n')
        assert 'valid=1' in run('bfs', small + ['--input', str(graph)])
        # Node 0 is isolated; the file's source 2 has a reachable neighbor.
        graph.write_text('3\n0 0\n0 0\n0 1\n2\n1\n1 1\n')
        assert 'valid=1' in run('bfs', small + ['--input', str(graph)])
        graph.write_text('3\n0 1\n1 0\n1 0\n0\n1\n9 1\n')
        run('bfs', small + ['--input', str(graph)], success=False)
        run('bfs', small + ['--input', str(empty)], success=False)
        assert 'valid=1' in run('pagerank', ['--size', '1', '--repeats', '1', '--warmup', '0'])
    print(f'PASS: {checks} CLI/input checks ({args.device})')


if __name__ == '__main__':
    main()
