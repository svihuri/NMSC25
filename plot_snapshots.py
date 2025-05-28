import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
plot_snapshots.py

This code reads the output of main.cpp simulation from run/output/ and creates a visualisation of the dye field at different time steps.
Output run/output/snapshot_*.png

Requires:
- pandas
- numpy
- matplotlib
"""

def load_data():
    """
    Load the data from the output file. Note hardcoded paths.
    return: pandas data frame
    """
    df = pd.read_csv('../run/output/output.csv', sep=';', usecols=[0,1,2,3,6])      # reading the whole file is too much at least for my laptop
    return df

def plot_snapshot(df, step, n, time=None):
    """
    Plot a snapshot of the dye field at a given time step.
    df: pandas data frame
    step: time step
    n: grid size
    """
    subset = df[df['step'] == step]
    D_grid = np.zeros((n, n))
    D_grid[subset['j'], subset['i']] = subset['D'].values

    fig, ax = plt.subplots()
    img = ax.imshow(np.zeros((n, n)), origin='lower', cmap='viridis', vmin=0, vmax=1)
    img.set_data(D_grid)
    ax.set_title(f"Time: {time} | Step: {step}")
    ax.set_xlabel('x')  
    ax.set_ylabel('y')

    plt.savefig(f'../run/output/snapshot_{step}.png')

def main():
    df = load_data()
    timesteps = [20, 100, 300, 500]     # timesteps to plot
    
    # grid size from data
    n = df['i'].max() + 1

    # pair step with time
    step_time_map = df[['step', 'time']].drop_duplicates().set_index('step')['time'].to_dict()

    for step in timesteps:
        time_val = step_time_map.get(step, None)
        plot_snapshot(df, step, n, time_val)

    print(f"Snapshots from simulation saved to ../run/output/snapshot_*.png")

if __name__ == "__main__":
    main()