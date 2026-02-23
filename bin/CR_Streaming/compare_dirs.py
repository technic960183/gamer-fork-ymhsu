#!/usr/bin/env python3
"""
Compare GAMER CR streaming outputs between two directories to verify
rotational symmetry: i.e., the same test problem run along X, Y, or Z
should produce identical results up to a permutation of vector components.

Usage:
    python3 compare_dirs.py [dir1] [dir2]

If no arguments are given, defaults to comparing
    bin/CR_Streaming  (auto-detected direction)
    bin/x_CR_Streaming (X-direction)
"""

import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────
#  Permutation tables
# ──────────────────────────────────────────────────────────────────────
# For a simulation whose streaming direction is along axis A, the
# "streaming-direction" component is labelled 1/X for A=X, 2/Y for A=Y,
# 3/Z for A=Z.  To compare with the X-direction reference we permute
# vector indices so that the streaming-dir component always maps to X.
#
# Convention:
#   X-sim  streaming=X  →  identity
#   Y-sim  streaming=Y  →  Y↔X, Z stays
#   Z-sim  streaming=Z  →  Z→X, X→Y, Y→Z

# Maps: src_direction -> { src_column : canonical_X_column }
VECTOR_PERMUTATION = {
    'X': {
        'MomX': 'MomX', 'MomY': 'MomY', 'MomZ': 'MomZ',
        'CR_F1': 'CR_F1', 'CR_F2': 'CR_F2', 'CR_F3': 'CR_F3',
        'ADV_VX': 'ADV_VX', 'ADV_VY': 'ADV_VY', 'ADV_VZ': 'ADV_VZ',
        'MagX': 'MagX', 'MagY': 'MagY', 'MagZ': 'MagZ',
    },
    'Y': {
        'MomY': 'MomX', 'MomX': 'MomY', 'MomZ': 'MomZ',
        'CR_F2': 'CR_F1', 'CR_F1': 'CR_F2', 'CR_F3': 'CR_F3',
        'ADV_VY': 'ADV_VX', 'ADV_VX': 'ADV_VY', 'ADV_VZ': 'ADV_VZ',
        'MagY': 'MagX', 'MagX': 'MagY', 'MagZ': 'MagZ',
    },
    'Z': {
        'MomZ': 'MomX', 'MomX': 'MomY', 'MomY': 'MomZ',
        'CR_F3': 'CR_F1', 'CR_F1': 'CR_F2', 'CR_F2': 'CR_F3',
        'ADV_VZ': 'ADV_VX', 'ADV_VX': 'ADV_VY', 'ADV_VY': 'ADV_VZ',
        'MagZ': 'MagX', 'MagX': 'MagY', 'MagY': 'MagZ',
    },
}

# Scalars: compared directly (no permutation)
SCALAR_VARS = ['Dens', 'Engy', 'CR_E', 'ADV_SIGMA', 'CRay', 'MagEngy']


# ──────────────────────────────────────────────────────────────────────
#  Data reader
# ──────────────────────────────────────────────────────────────────────
class GAMERLineData:
    """Parse a GAMER Xline / Yline / Zline output file."""

    def __init__(self, filename, record_dump_file=None):
        self.filename = filename
        self.time = None
        self.step = None
        self.columns = []
        self.data = self._read_file()
        if record_dump_file and os.path.exists(record_dump_file):
            self._read_time_info(record_dump_file)

    def _read_file(self):
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        header = lines[0]
        self.columns = [c for c in header.strip().split() if c != '#']

        data = {col: [] for col in self.columns}
        for line in lines[1:]:
            if line.strip() and not line.startswith('#'):
                vals = line.strip().split()
                for i, col in enumerate(self.columns):
                    if i < len(vals):
                        try:
                            data[col].append(float(vals[i]))
                        except ValueError:
                            pass
        for col in data:
            data[col] = np.array(data[col])
        return data

    def _read_time_info(self, record_dump_file):
        dump_id = int(os.path.basename(self.filename).split('_')[-1])
        try:
            with open(record_dump_file, 'r') as f:
                for line in f.readlines()[1:]:
                    parts = line.strip().split()
                    if parts and int(parts[0]) == dump_id:
                        self.time = float(parts[1])
                        self.step = int(parts[2])
                        break
        except (IOError, ValueError):
            pass

    def get(self, varname):
        return self.data.get(varname, None)


# ──────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────
def detect_direction(directory):
    """Return ('X','Y', or 'Z'), the line-type prefix, and the
    coordinate column name for files found in *directory*."""
    for prefix, coord in [('Xline', 'x'), ('Yline', 'y'), ('Zline', 'z')]:
        if glob.glob(os.path.join(directory, '{}_*_000000'.format(prefix))):
            return prefix[0], prefix, coord
    return None, None, None


def get_sorted_files(directory, prefix):
    """Return a list of (filename, timestep_str) sorted by timestep."""
    pattern = os.path.join(directory, '{}_*'.format(prefix))
    files = sorted(glob.glob(pattern))
    result = []
    for f in files:
        ts = os.path.basename(f).split('_')[-1]
        result.append((f, ts))
    return result


# ──────────────────────────────────────────────────────────────────────
#  Comparator
# ──────────────────────────────────────────────────────────────────────
class SymmetryComparator:
    """Compare two GAMER line-output directories with proper
    permutation of vector components."""

    def __init__(self, dir1, dir2, label1=None, label2=None):
        self.dir1 = dir1
        self.dir2 = dir2

        # Detect directions
        self.dir1_axis, self.dir1_prefix, self.dir1_coord = detect_direction(dir1)
        self.dir2_axis, self.dir2_prefix, self.dir2_coord = detect_direction(dir2)

        if self.dir1_axis is None:
            raise RuntimeError("No Xline/Yline/Zline files in " + dir1)
        if self.dir2_axis is None:
            raise RuntimeError("No Xline/Yline/Zline files in " + dir2)

        self.label1 = label1 or os.path.basename(dir1)
        self.label2 = label2 or os.path.basename(dir2)

        # Get permutation tables
        self.perm1 = VECTOR_PERMUTATION[self.dir1_axis]
        self.perm2 = VECTOR_PERMUTATION[self.dir2_axis]

        # Build comparison pairs: (canonical_name, dir1_col, dir2_col)
        # The canonical name is the X-direction label.
        # perm maps  src_col -> canonical(X) col,
        # so inverse maps  canonical(X) col -> src_col.
        inv_perm1 = {v: k for k, v in self.perm1.items()}
        inv_perm2 = {v: k for k, v in self.perm2.items()}

        self.compare_pairs = []

        # Scalars
        for s in SCALAR_VARS:
            self.compare_pairs.append((s, s, s))

        # Vectors
        for canon in ['MomX', 'MomY', 'MomZ',
                       'CR_F1', 'CR_F2', 'CR_F3',
                       'ADV_VX', 'ADV_VY', 'ADV_VZ',
                       'MagX', 'MagY', 'MagZ']:
            col1 = inv_perm1.get(canon, canon)
            col2 = inv_perm2.get(canon, canon)
            self.compare_pairs.append((canon, col1, col2))

    # ── file pairing ────────────────────────────────────────────────
    def get_file_pairs(self):
        files1 = get_sorted_files(self.dir1, self.dir1_prefix)
        files2 = get_sorted_files(self.dir2, self.dir2_prefix)
        ts_map2 = {ts: f for f, ts in files2}
        pairs = []
        for f1, ts in files1:
            if ts in ts_map2:
                pairs.append((f1, ts_map2[ts], ts))
        return pairs

    # ── single-timestep comparison ──────────────────────────────────
    def compare_timestep(self, file1, file2, timestep):
        rd1 = os.path.join(self.dir1, 'Record__Dump')
        rd2 = os.path.join(self.dir2, 'Record__Dump')
        d1 = GAMERLineData(file1, rd1)
        d2 = GAMERLineData(file2, rd2)

        print("\n" + "=" * 80)
        print("Timestep {}".format(timestep))
        print("=" * 80)
        print("  {} ({}-dir): {}".format(self.label1, self.dir1_axis,
                                          os.path.basename(file1)))
        print("  {} ({}-dir): {}".format(self.label2, self.dir2_axis,
                                          os.path.basename(file2)))
        if d1.time is not None:
            print("  {} : t = {:.6e}, step = {}".format(self.label1, d1.time, d1.step))
        if d2.time is not None:
            print("  {} : t = {:.6e}, step = {}".format(self.label2, d2.time, d2.step))
        print("-" * 80)

        # Coordinates along the streaming direction
        coord1 = d1.get(self.dir1_coord)   # e.g. 'y' for Y-sim
        coord2 = d2.get(self.dir2_coord)   # e.g. 'x' for X-sim

        results = {
            'time1': d1.time, 'step1': d1.step,
            'time2': d2.time, 'step2': d2.step,
            'coord1': coord1, 'coord2': coord2,
        }

        # Check coordinate alignment
        if coord1 is not None and coord2 is not None:
            if len(coord1) == len(coord2) and np.allclose(coord1, coord2, rtol=1e-10):
                print("  Coordinates ({} vs {}): match".format(
                    self.dir1_coord, self.dir2_coord))
            else:
                print("  Coordinates ({} vs {}): *** MISMATCH ***  "
                      "len {} vs {}, range [{:.6f},{:.6f}] vs [{:.6f},{:.6f}]".format(
                          self.dir1_coord, self.dir2_coord,
                          len(coord1), len(coord2),
                          coord1[0], coord1[-1], coord2[0], coord2[-1]))

        print()
        print("  {:<20s} {:<20s} {:<20s} {:>12s} {:>12s} {:>12s}".format(
            "Canonical", self.label1 + " col", self.label2 + " col",
            "Max|Diff|", "Mean|Diff|", "MaxRelDiff"))
        print("  " + "-" * 96)

        for canon, col1, col2 in self.compare_pairs:
            v1 = d1.get(col1)
            v2 = d2.get(col2)

            if v1 is None or v2 is None:
                tag = ""
                if v1 is None:
                    tag += " (missing in {})".format(self.label1)
                if v2 is None:
                    tag += " (missing in {})".format(self.label2)
                print("  {:<20s} {:<20s} {:<20s}  SKIP{}".format(
                    canon, col1, col2, tag))
                continue

            n = min(len(v1), len(v2))
            v1, v2 = v1[:n], v2[:n]

            diff = v1 - v2
            abs_diff = np.abs(diff)
            denom = np.maximum(np.abs(v1), np.abs(v2))
            rel_diff = np.where(denom > 0, abs_diff / denom, 0.0)

            max_abs = np.max(abs_diff)
            mean_abs = np.mean(abs_diff)
            max_rel = np.max(rel_diff)

            flag = ""
            if max_abs == 0:
                flag = "  identical"
            elif max_rel < 1e-12:
                flag = "  ~machine-eps"
            elif max_rel > 0.01:
                flag = "  *** LARGE ***"

            print("  {:<20s} {:<20s} {:<20s} {:>12.4e} {:>12.4e} {:>12.4e}{}".format(
                canon, col1, col2, max_abs, mean_abs, max_rel, flag))

            results[canon] = {
                'col1': col1, 'col2': col2,
                'data1': v1, 'data2': v2,
                'diff': diff, 'abs_diff': abs_diff, 'rel_diff': rel_diff,
                'max_abs': max_abs, 'mean_abs': mean_abs, 'max_rel': max_rel,
            }

        return results

    # ── run all timesteps ───────────────────────────────────────────
    def compare_all(self, save_plots=True):
        pairs = self.get_file_pairs()
        if not pairs:
            print("No matching timestep pairs found!")
            return None

        print("\n{} ({}-dir) vs {} ({}-dir)".format(
            self.label1, self.dir1_axis, self.label2, self.dir2_axis))
        print("Found {} matching timesteps".format(len(pairs)))
        print("\nVector component mapping (after permutation):")
        for canon, col1, col2 in self.compare_pairs:
            if col1 != col2:
                print("  {} col '{}' <-> {} col '{}'  (canonical: {})".format(
                    self.label1, col1, self.label2, col2, canon))

        all_results = {}
        for f1, f2, ts in pairs:
            all_results[ts] = self.compare_timestep(f1, f2, ts)

        if save_plots:
            self._plot_profiles(all_results, pairs)
            self._plot_diff_evolution(all_results, pairs)

        return all_results

    # ── profile plots ───────────────────────────────────────────────
    def _plot_profiles(self, all_results, pairs):
        print("\n" + "=" * 80)
        print("Generating comparison plots...")
        print("=" * 80)

        plot_vars = [c for c in ['Dens', 'CR_E', 'CR_F1', 'CR_F2', 'CR_F3', 'ADV_SIGMA',
                                  'MagX', 'MomX']
                     if any(c in all_results[ts] for _, _, ts in pairs)]

        if not plot_vars:
            return

        for var in plot_vars:
            n_ts = min(len(pairs), 4)
            fig, axes = plt.subplots(n_ts, 2, figsize=(14, 4 * n_ts),
                                      squeeze=False)
            # Find the actual column names used
            sample_ts = pairs[0][2]
            if var in all_results[sample_ts]:
                col1_name = all_results[sample_ts][var]['col1']
                col2_name = all_results[sample_ts][var]['col2']
            else:
                col1_name = col2_name = var

            fig.suptitle("Symmetry: canonical '{}'\n{} col '{}' vs {} col '{}'".format(
                var, self.label1, col1_name, self.label2, col2_name),
                fontsize=13, y=1.02)

            for idx in range(n_ts):
                f1, f2, ts = pairs[idx]
                res = all_results[ts]

                if var not in res:
                    continue

                col1 = res[var]['col1']
                col2 = res[var]['col2']
                v1 = res[var]['data1']
                v2 = res[var]['data2']

                coord1 = res['coord1']
                coord2 = res['coord2']

                # Left panel: overlay profiles
                ax = axes[idx, 0]
                if coord1 is not None and coord2 is not None:
                    n = min(len(coord1), len(v1))
                    ax.plot(coord1[:n], v1[:n], 'b-', lw=2,
                            label='{} [{}]'.format(self.label1, col1))
                    n2 = min(len(coord2), len(v2))
                    ax.plot(coord2[:n2], v2[:n2], 'r--', lw=2,
                            label='{} [{}]'.format(self.label2, col2))
                else:
                    ax.plot(v1, 'b-', lw=2, label=self.label1)
                    ax.plot(v2, 'r--', lw=2, label=self.label2)

                title = 'ts={}'.format(ts)
                t = res.get('time1')
                if t is not None:
                    title += ', t={:.3e}'.format(t)
                ax.set_title(title)
                ax.set_ylabel(var)
                ax.set_xlabel('Position along streaming dir')
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)

                # Right panel: difference
                ax2 = axes[idx, 1]
                x = coord1 if coord1 is not None else np.arange(len(v1))
                n = min(len(x), len(res[var]['abs_diff']))
                ax2.semilogy(x[:n], res[var]['abs_diff'][:n] + 1e-30,
                             'm-', lw=1.5, label='|diff|')
                ax2.semilogy(x[:n], res[var]['rel_diff'][:n] + 1e-30,
                             'g--', lw=1.5, label='rel diff')
                ax2.set_title('Difference (ts={})'.format(ts))
                ax2.set_ylabel('Difference')
                ax2.set_xlabel('Position along streaming dir')
                ax2.legend(fontsize=8)
                ax2.grid(True, alpha=0.3)

            plt.tight_layout()
            outf = os.path.join(self.dir1, 'symm_cmp_{}.png'.format(var))
            plt.savefig(outf, dpi=150, bbox_inches='tight')
            print("  Saved: {}".format(outf))
            plt.close()

    # ── difference-vs-time plot ─────────────────────────────────────
    def _plot_diff_evolution(self, all_results, pairs):
        plot_vars = [c for c in ['Dens', 'CR_E', 'CR_F1', 'ADV_SIGMA',
                                  'MagX', 'MomX']
                     if any(c in all_results[ts] for _, _, ts in pairs)]
        if not plot_vars:
            return

        fig, axes = plt.subplots(len(plot_vars), 2,
                                  figsize=(12, 3.5 * len(plot_vars)),
                                  squeeze=False)
        fig.suptitle('{} ({}) vs {} ({})  -- error evolution'.format(
            self.label1, self.dir1_axis, self.label2, self.dir2_axis),
            fontsize=14, y=1.01)

        for vi, var in enumerate(plot_vars):
            ts_list, max_abs, mean_abs, max_rel = [], [], [], []
            for _, _, ts in pairs:
                if var in all_results[ts]:
                    ts_list.append(int(ts))
                    max_abs.append(all_results[ts][var]['max_abs'])
                    mean_abs.append(all_results[ts][var]['mean_abs'])
                    max_rel.append(all_results[ts][var]['max_rel'])

            if not ts_list:
                continue

            ax = axes[vi, 0]
            ax.semilogy(ts_list, max_abs, 'ro-', label='max |diff|', ms=6)
            ax.semilogy(ts_list, mean_abs, 'bs-', label='mean |diff|', ms=6)
            ax.set_ylabel('Absolute difference')
            ax.set_xlabel('Dump ID')
            ax.set_title(var)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

            ax = axes[vi, 1]
            ax.semilogy(ts_list, max_rel, 'ro-', label='max rel diff', ms=6)
            ax.set_ylabel('Relative difference')
            ax.set_xlabel('Dump ID')
            ax.set_title(var)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        outf = os.path.join(self.dir1, 'symm_cmp_evolution.png')
        plt.savefig(outf, dpi=150, bbox_inches='tight')
        print("  Saved: {}".format(outf))
        plt.close()


# ──────────────────────────────────────────────────────────────────────
#  Main
# ──────────────────────────────────────────────────────────────────────
def main():
    base = '/global/u1/y/ymhsu/gamer-fork-ymhsu/bin'

    if len(sys.argv) >= 3:
        dir1 = sys.argv[1]
        dir2 = sys.argv[2]
    else:
        dir1 = os.path.join(base, 'CR_Streaming')
        dir2 = os.path.join(base, 'x_CR_Streaming')

    for d in [dir1, dir2]:
        if not os.path.isdir(d):
            print("ERROR: directory not found: {}".format(d))
            sys.exit(1)

    label1 = os.path.basename(dir1)
    label2 = os.path.basename(dir2)

    comp = SymmetryComparator(dir1, dir2, label1, label2)
    results = comp.compare_all(save_plots=True)

    print("\n" + "=" * 80)
    print("Done!  Plots saved to {}".format(dir1))
    print("=" * 80)


if __name__ == '__main__':
    main()
