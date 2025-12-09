#!/usr/bin/env python3
"""
Visualize multipole simulation output
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from pathlib import Path


def visualize_frames(output_dir='output', save_animation=False):
    """
    Visualize simulation frames
    """

    # Find all frame files
    frame_files = sorted(glob.glob(f"{output_dir}/frame_*.csv"))

    if len(frame_files) == 0:
        print(f"No frames found in {output_dir}/")
        return

    print(f"Found {len(frame_files)} frames")

    # Read first frame to get bounds
    data = np.genfromtxt(frame_files[0], delimiter=',', skip_header=1)
    x_all = data[:, 0]
    y_all = data[:, 1]

    # Set up plot
    plt.figure(figsize=(8, 8))

    # Animation or static frames
    if save_animation:
        print("Saving animation...")
        frames = []

    for i, frame_file in enumerate(frame_files):
        data = np.genfromtxt(frame_file, delimiter=',', skip_header=1)
        x = data[:, 0]
        y = data[:, 1]
        q = data[:, 2]

        plt.clf()

        # Color by charge
        plt.scatter(x, y, c=q, s=10, cmap='RdBu_r', vmin=-1, vmax=1, alpha=0.7)
        plt.colorbar(label='Charge')

        plt.xlim(-1.1, 1.1)
        plt.ylim(-1.1, 1.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'Frame {i} / {len(frame_files)-1}')
        plt.grid(True, alpha=0.3)
        plt.axis('equal')

        if save_animation:
            # Save frame for animation
            plt.savefig(f'{output_dir}/plot_{i:06d}.png', dpi=100)
        else:
            plt.pause(0.05)

    if save_animation:
        print(f"Frames saved to {output_dir}/plot_*.png")
        print("Create animation with:")
        print(f"  ffmpeg -r 20 -i {output_dir}/plot_%06d.png -c:v libx264 -vf fps=20 -pix_fmt yuv420p animation.mp4")
    else:
        plt.show()


def plot_energy_evolution(output_dir='output'):
    """
    Plot energy evolution over time (if energy data available)
    """
    # This would require saving energy data - placeholder for future
    pass


def analyze_distribution(output_dir='output'):
    """
    Analyze particle distribution
    """
    frame_files = sorted(glob.glob(f"{output_dir}/frame_*.csv"))

    if len(frame_files) == 0:
        print(f"No frames found in {output_dir}/")
        return

    # Compare initial and final distributions
    data_init = np.genfromtxt(frame_files[0], delimiter=',', skip_header=1)
    data_final = np.genfromtxt(frame_files[-1], delimiter=',', skip_header=1)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Initial
    axes[0].scatter(data_init[:, 0], data_init[:, 1], c=data_init[:, 2],
                   s=10, cmap='RdBu_r', vmin=-1, vmax=1, alpha=0.7)
    axes[0].set_xlim(-1.1, 1.1)
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_title('Initial Distribution')
    axes[0].grid(True, alpha=0.3)
    axes[0].axis('equal')

    # Final
    axes[1].scatter(data_final[:, 0], data_final[:, 1], c=data_final[:, 2],
                   s=10, cmap='RdBu_r', vmin=-1, vmax=1, alpha=0.7)
    axes[1].set_xlim(-1.1, 1.1)
    axes[1].set_ylim(-1.1, 1.1)
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    axes[1].set_title('Final Distribution')
    axes[1].grid(True, alpha=0.3)
    axes[1].axis('equal')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/distribution_comparison.png', dpi=150)
    print(f"Saved: {output_dir}/distribution_comparison.png")
    plt.show()


def main():
    """Main function"""

    import argparse

    parser = argparse.ArgumentParser(description='Visualize multipole simulation')
    parser.add_argument('--dir', default='output', help='Output directory')
    parser.add_argument('--save', action='store_true', help='Save animation frames')
    parser.add_argument('--analyze', action='store_true', help='Analyze distributions')

    args = parser.parse_args()

    if args.analyze:
        analyze_distribution(args.dir)
    else:
        visualize_frames(args.dir, args.save)


if __name__ == '__main__':
    main()
