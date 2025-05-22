import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

df = pd.read_csv("output.csv", sep=';') 

# grid size from data
n = df['i'].max() + 1

# timestep and time
steps = df['step'].unique()
steps.sort()
time = df['time'].unique()
time.sort()

fig, ax = plt.subplots()
img = ax.imshow(np.zeros((n, n)), origin='lower', cmap='viridis', vmin=0, vmax=1)

# update function for animation
def update(frame):
    step = steps[frame]
    subset = df[df['step'] == step]
    D_grid = np.zeros((n, n))
    D_grid[subset['i'], subset['j']] = subset['D'].values
    img.set_data(D_grid)
    ax.set_title(f"Time: {time[frame]} | Step: {step}")
    return [img]

# create the animation
ani = FuncAnimation(fig, update, frames=len(steps), interval=100, blit=True)
ani.save("D_animation.mp4", fps=10)






