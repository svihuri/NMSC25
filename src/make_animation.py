import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""
make_animation.py
This code reads the output of main.cpp simulation from run/output/ and creates an animation of the dye field.
Output run/output/simulation.gif

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


def main():
    df = load_data()
    # grid size from data
    n = df['i'].max() + 1

    # timestep and time
    steps = df['step'].unique()
    steps.sort()
    time = df['time'].unique()
    time.sort()

    fig, ax = plt.subplots()
    img = ax.imshow(np.zeros((n, n)), origin='lower', cmap='viridis', vmin=0, vmax=1)

    def update(frame):
        """
        Update the image for the animation.
        frame: frame number
        return: AxesImage
        """
        step = steps[frame]
        subset = df[df['step'] == step]
        D_grid = np.zeros((n, n))
        D_grid[subset['j'], subset['i']] = subset['D'].values
        img.set_data(D_grid)
        ax.set_title(f"Time: {time[frame]} | Step: {step}")
        return [img]

    # create the animation
    ani = FuncAnimation(fig, update, frames=len(steps), interval=100, blit=True)
    ani.save("../run/output/simulation.gif", fps=12)

    print("Animation saved to ../run/output/simulation.gif")

if __name__ == "__main__":
    main()






