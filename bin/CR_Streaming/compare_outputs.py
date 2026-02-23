#!/usr/bin/env python3
"""
Compare GAMER and Athena++ CR streaming outputs step by step.
This script reads both output formats and compares key variables to identify differences.
"""

import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from pathlib import Path


class GAMERData:
    """Parse GAMER Xline output files"""

    def __init__(self, filename, record_dump_file=None):
        self.filename = filename
        self.time = None
        self.step = None
        self.data = self._read_file()
        if record_dump_file:
            self._read_time_info(record_dump_file)

    def _read_file(self):
        """Read GAMER output file and parse into structured data"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        # Parse header to get column names
        header_line = lines[0]
        columns = header_line.strip().split()
        columns = [col for col in columns if col != '#']

        # Read data
        data_dict = {col: [] for col in columns}

        for line in lines[1:]:
            if line.strip() and not line.startswith('#'):
                values = line.strip().split()
                for i, col in enumerate(columns):
                    if i < len(values):
                        try:
                            data_dict[col].append(float(values[i]))
                        except ValueError:
                            continue

        # Convert to numpy arrays
        for col in data_dict:
            data_dict[col] = np.array(data_dict[col])

        return data_dict

    def _read_time_info(self, record_dump_file):
        """Read time and step from Record__Dump file"""
        # Extract DumpID from filename (e.g., Xline_y0.000_z0.000_000000 -> 0)
        dump_id = int(self.filename.split('_')[-1])
        
        try:
            with open(record_dump_file, 'r') as f:
                lines = f.readlines()
            
            # Skip header line
            for line in lines[1:]:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 3 and int(parts[0]) == dump_id:
                        self.time = float(parts[1])
                        self.step = int(parts[2])
                        break
        except (IOError, ValueError) as e:
            print("Warning: Could not read time info from {}: {}".format(record_dump_file, e))

    def get_variable(self, varname):
        """Get a specific variable"""
        return self.data.get(varname, None)


class AthenaData:
    """Parse Athena++ tabular output files"""

    def __init__(self, filename):
        self.filename = filename
        self.time = None
        self.cycle = None
        self.data = self._read_file()

    def _read_file(self):
        """Read Athena++ output file and parse into structured data"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()

        # Parse first line for time and cycle info
        # Format: # Athena++ data at time=0.000000e+00  cycle=0  variables=prim
        if lines and lines[0].startswith('#'):
            first_line = lines[0]
            if 'time=' in first_line and 'cycle=' in first_line:
                try:
                    time_part = first_line.split('time=')[1].split()[0]
                    cycle_part = first_line.split('cycle=')[1].split()[0]
                    self.time = float(time_part)
                    self.cycle = int(cycle_part)
                except (ValueError, IndexError) as e:
                    print("Warning: Could not parse time/cycle from {}: {}".format(self.filename, e))

        # Find header line (starts with # and contains variable names)
        header_line = None
        data_start_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('#') and 'x1v' in line:
                header_line = line
                data_start_idx = i + 1
                break

        if header_line is None:
            raise ValueError("Could not find header in " + self.filename)

        # Parse column names
        columns = header_line.strip().replace('#', '').split()

        # Read data
        data_dict = {col: [] for col in columns}

        for line in lines[data_start_idx:]:
            if line.strip() and not line.startswith('#'):
                values = line.strip().split()
                # First column is 'i' (index), values[0] -> columns[0]
                for i, col in enumerate(columns):
                    if i < len(values):
                        try:
                            data_dict[col].append(float(values[i]))
                        except (ValueError, IndexError):
                            continue

        # Convert to numpy arrays
        for col in data_dict:
            data_dict[col] = np.array(data_dict[col])

        return data_dict

    def get_variable(self, varname):
        """Get a specific variable"""
        return self.data.get(varname, None)


class OutputComparator:
    """Compare GAMER and Athena outputs"""

    def __init__(self, gamer_dir, athena_dir):
        self.gamer_dir = gamer_dir
        self.athena_dir = athena_dir

        # Variable mappings between GAMER and Athena
        self.var_mapping = {
            'x': 'x1v',           # Position
            'Dens': 'rho',        # Density
            'Engy': None,         # Total energy (computed differently)
            'CR_E': 'Ec',         # CR energy density
            'CR_F1': 'Fc1',       # CR flux in x
            'CR_F2': 'Fc2',       # CR flux in y
            'CR_F3': 'Fc3',       # CR flux in z
            'ADV_SIGMA': 'Sigma_adv1',  # Advection sigma
            'ADV_VX': 'Vc1',      # Advection velocity x
            'ADV_VY': 'Vc2',      # Advection velocity y
            'ADV_VZ': 'Vc3',      # Advection velocity z
            'MagX': 'Bcc1',       # Magnetic field x
            'MagY': 'Bcc2',       # Magnetic field y
            'MagZ': 'Bcc3',       # Magnetic field z
        }

    def get_file_pairs(self):
        """Get matching pairs of GAMER and Athena output files"""
        gamer_files = sorted(glob.glob(os.path.join(self.gamer_dir, 'Xline_y0.000_z0.000_0*')))
        athena_files = sorted(glob.glob(os.path.join(self.athena_dir, 'cr.block0.out1.0*.tab')))

        # Extract timestep numbers - GAMER has 6 digits, Athena has 5
        pairs = []
        for gf in gamer_files:
            gamer_timestep = gf.split('_')[-1]  # e.g., "000000"
            # Convert to Athena format (5 digits)
            athena_timestep = gamer_timestep[-5:]  # Last 5 digits
            # Find matching Athena file
            af_pattern = os.path.join(self.athena_dir, 'cr.block0.out1.' + athena_timestep + '.tab')
            if os.path.exists(af_pattern):
                pairs.append((gf, af_pattern, gamer_timestep))

        return pairs

    def compare_timestep(self, gamer_file, athena_file, timestep):
        """Compare a single timestep between GAMER and Athena"""
        print("\n" + "="*80)
        print("Comparing Timestep {}".format(timestep))
        print("="*80)
        print("GAMER:  {}".format(os.path.basename(gamer_file)))
        print("Athena: {}".format(os.path.basename(athena_file)))
        print("-"*80)

        # Load data
        record_dump = os.path.join(self.gamer_dir, 'Record__Dump')
        gamer = GAMERData(gamer_file, record_dump if os.path.exists(record_dump) else None)
        athena = AthenaData(athena_file)
        
        # Print time/step information
        if gamer.time is not None and gamer.step is not None:
            print("GAMER:  Time = {:.6e}, Step = {}".format(gamer.time, gamer.step))
        if athena.time is not None and athena.cycle is not None:
            print("Athena: Time = {:.6e}, Cycle = {}".format(athena.time, athena.cycle))
        print("-"*80)

        # Compare each variable
        results = {
            'gamer_time': gamer.time,
            'gamer_step': gamer.step,
            'athena_time': athena.time,
            'athena_cycle': athena.cycle,
        }

        for gamer_var, athena_var in self.var_mapping.items():
            if athena_var is None:
                continue

            g_data = gamer.get_variable(gamer_var)
            a_data = athena.get_variable(athena_var)

            if g_data is None or a_data is None:
                print("  {:15s} -> {:15s}: MISSING DATA".format(gamer_var, athena_var))
                continue

            # Ensure same size
            min_len = min(len(g_data), len(a_data))
            if min_len == 0:
                print("  {:15s} -> {:15s}: EMPTY DATA".format(gamer_var, athena_var))
                continue

            g_data = g_data[:min_len]
            a_data = a_data[:min_len]

            # If x coordinates, adjust for box offset
            if gamer_var == 'x':
                g_data = g_data - 1  # Adjust for box offset

            # Compute statistics
            diff = g_data - a_data
            abs_diff = np.abs(diff)
            rel_diff = abs_diff / (np.abs(a_data) + 1e-30)

            p25_rel_diff = np.percentile(rel_diff, 25)
            p50_rel_diff = np.percentile(rel_diff, 50)
            p75_rel_diff = np.percentile(rel_diff, 75)
            
            p25_abs_diff = np.percentile(abs_diff, 25)
            p50_abs_diff = np.percentile(abs_diff, 50)
            p75_abs_diff = np.percentile(abs_diff, 75)

            results[gamer_var] = {
                'gamer': g_data,
                'athena': a_data,
                'diff': diff,
                'abs_diff': abs_diff,
                'rel_diff': rel_diff,
                'p25_rel_diff': p25_rel_diff,
                'p50_rel_diff': p50_rel_diff,
                'p75_rel_diff': p75_rel_diff,
                'p25_abs_diff': p25_abs_diff,
                'p50_abs_diff': p50_abs_diff,
                'p75_abs_diff': p75_abs_diff,
            }

            # Print comparison
            print("  {:15s} -> {:15s}:".format(gamer_var, athena_var))
            print("    P25 Rel Diff: {:.6e}  P50 Rel Diff: {:.6e}  P75 Rel Diff: {:.6e}".format(p25_rel_diff, p50_rel_diff, p75_rel_diff))

            # Flag significant differences
            if p75_rel_diff > 0.01:  # 1% difference
                print("    *** WARNING: Large relative difference! ***")

        return results

    def compare_all(self, save_plots=True):
        """Compare all timesteps and optionally generate plots"""
        pairs = self.get_file_pairs()

        if not pairs:
            print("No matching file pairs found!")
            return

        print("\nFound {} matching timestep pairs".format(len(pairs)))

        all_results = {}

        for gamer_file, athena_file, timestep in pairs:
            results = self.compare_timestep(gamer_file, athena_file, timestep)
            all_results[timestep] = results

        if save_plots:
            self.generate_comparison_plots(all_results, pairs)

        return all_results

    def generate_comparison_plots(self, all_results, pairs):
        """Generate comparison plots for key variables"""
        print("\n" + "="*80)
        print("Generating comparison plots...")
        print("="*80)

        # Variables to plot
        key_vars = ['CR_E', 'CR_F1', 'ADV_SIGMA']

        for var in key_vars:
            if var not in self.var_mapping or self.var_mapping[var] is None:
                continue

            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle('Comparison: {} vs {}'.format(var, self.var_mapping[var]), fontsize=16)

            for idx, (gamer_file, athena_file, timestep) in enumerate(pairs):
                if idx >= 4:  # Only plot first 4 timesteps
                    break

                ax = axes[idx // 2, idx % 2]

                results = all_results[timestep]
                if var not in results:
                    continue

                x_gamer = all_results[timestep].get('x', {}).get('gamer', None)
                x_athena = all_results[timestep].get('x', {}).get('athena', None)

                # if x_gamer is not None:
                #     x_gamer = x_gamer - 1  # Adjust for box offset

                g_data = results[var]['gamer']
                a_data = results[var]['athena']
                abs_diff = results[var]['abs_diff']
                rel_diff = results[var]['rel_diff']

                # Check if x coordinates match
                x_match = True
                if x_gamer is not None and x_athena is not None:
                    if len(x_gamer) != len(x_athena) or not np.allclose(x_gamer, x_athena, rtol=1e-5):
                        print("  *** WARNING: x-coordinates do not match for {} at timestep {}!".format(var, timestep))
                        print("      GAMER x range: [{:.6f}, {:.6f}], length: {}".format(
                            x_gamer[0], x_gamer[-1], len(x_gamer)))
                        print("      Athena x range: [{:.6f}, {:.6f}], length: {}".format(
                            x_athena[0], x_athena[-1], len(x_athena)))
                        x_match = False

                # Plot data on primary axis
                if x_gamer is not None and x_athena is not None:
                    ax.plot(x_gamer, g_data, 'b-', label='GAMER', linewidth=2)
                    ax.plot(x_athena, a_data, 'r--', label='Athena', linewidth=2)
                    x_axis = x_gamer
                else:
                    ax.plot(g_data, 'b-', label='GAMER', linewidth=2)
                    ax.plot(a_data, 'r--', label='Athena', linewidth=2)
                    x_axis = np.arange(len(g_data))

                # Create title with time/step info
                title_parts = ['Timestep {}'.format(timestep)]
                gamer_time = results.get('gamer_time')
                gamer_step = results.get('gamer_step')
                if gamer_time is not None:
                    title_parts.append('t={:.3e}'.format(gamer_time))
                if gamer_step is not None:
                    title_parts.append('step={}'.format(gamer_step))
                ax.set_title(', '.join(title_parts))
                ax.set_xlabel('Position' if x_gamer is not None else 'Index')
                ax.set_ylabel(var, color='b')
                ax.tick_params(axis='y')
                ax.legend(loc='upper left')
                ax.grid(True, alpha=0.3)

                # Create secondary axis for errors only if x coordinates match
                if x_match:
                    ax2 = ax.twinx()
                    ax2.plot(x_axis, rel_diff, 'g-', label='Relative Error', linewidth=1.5, alpha=0.7)
                    ax2.plot(x_axis, abs_diff, 'm--', label='Absolute Error', linewidth=1.5, alpha=0.7)
                    ax2.set_ylabel('Error (log scale)', color='g')
                    ax2.set_yscale('log')
                    ax2.tick_params(axis='y')
                    ax2.legend(loc='upper right')

            plt.tight_layout()
            output_file = os.path.join(self.gamer_dir, 'comparison_{}.png'.format(var))
            plt.savefig(output_file, dpi=150)
            print("  Saved: {}".format(output_file))
            plt.close()

        # Also create a difference plot
        self._plot_differences(all_results, pairs)

    def _plot_differences(self, all_results, pairs):
        """Plot differences over time for key variables"""
        key_vars = ['CR_E', 'CR_F1']

        fig, axes = plt.subplots(len(key_vars), 2, figsize=(12, 4*len(key_vars)))
        if len(key_vars) == 1:
            axes = axes.reshape(1, -1)

        for idx, var in enumerate(key_vars):
            ax_rel = axes[idx, 0]
            ax_abs = axes[idx, 1]

            timesteps = []
            p25_rel_diffs = []
            p50_rel_diffs = []
            p75_rel_diffs = []
            p25_abs_diffs = []
            p50_abs_diffs = []
            p75_abs_diffs = []

            for gamer_file, athena_file, timestep in pairs:
                results = all_results[timestep]
                if var in results:
                    timesteps.append(int(timestep))
                    p25_rel_diffs.append(results[var]['p25_rel_diff'])
                    p50_rel_diffs.append(results[var]['p50_rel_diff'])
                    p75_rel_diffs.append(results[var]['p75_rel_diff'])
                    p25_abs_diffs.append(results[var]['p25_abs_diff'])
                    p50_abs_diffs.append(results[var]['p50_abs_diff'])
                    p75_abs_diffs.append(results[var]['p75_abs_diff'])

            if timesteps:
                # Plot relative differences
                ax_rel.plot(timesteps, p25_rel_diffs, 'go-', label='P25', linewidth=2, markersize=8)
                ax_rel.plot(timesteps, p50_rel_diffs, 'bs-', label='P50', linewidth=2, markersize=8)
                ax_rel.plot(timesteps, p75_rel_diffs, 'ro-', label='P75', linewidth=2, markersize=8)
                ax_rel.set_ylabel('Relative Difference')
                ax_rel.set_xlabel('Timestep')
                ax_rel.set_title('{} - Relative Differences over Time'.format(var))
                ax_rel.legend()
                ax_rel.grid(True, alpha=0.3)
                ax_rel.set_yscale('log')
                
                # Plot absolute differences
                ax_abs.plot(timesteps, p25_abs_diffs, 'go-', label='P25', linewidth=2, markersize=8)
                ax_abs.plot(timesteps, p50_abs_diffs, 'bs-', label='P50', linewidth=2, markersize=8)
                ax_abs.plot(timesteps, p75_abs_diffs, 'ro-', label='P75', linewidth=2, markersize=8)
                ax_abs.set_ylabel('Absolute Difference')
                ax_abs.set_xlabel('Timestep')
                ax_abs.set_title('{} - Absolute Differences over Time'.format(var))
                ax_abs.legend()
                ax_abs.grid(True, alpha=0.3)
                ax_abs.set_yscale('log')

        plt.tight_layout()
        output_file = os.path.join(self.gamer_dir, 'comparison_differences_vs_time.png')
        plt.savefig(output_file, dpi=150)
        print("  Saved: {}".format(output_file))
        plt.close()


def main():
    """Main function"""
    # Set directories
    gamer_dir = '/global/u1/y/ymhsu/gamer-fork-ymhsu/bin/CR_Streaming'
    athena_dir = '/global/u1/y/ymhsu/athena/working'

    # Create comparator
    comparator = OutputComparator(gamer_dir, athena_dir)

    # Run comparison
    results = comparator.compare_all(save_plots=True)

    print("\n" + "="*80)
    print("Comparison Complete!")
    print("="*80)
    print("\nResults saved to: {}".format(gamer_dir))
    print("  - comparison_*.png: Variable comparisons for each timestep")
    print("  - comparison_differences_vs_time.png: Difference evolution")


if __name__ == '__main__':
    main()
