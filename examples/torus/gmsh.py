import sys, os, re
import copy
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
import pymech as pm
from gmsh_writer import generate_gmsh_script

if __name__ == "__main__":
   # Given parameters
   half = False
   R = 1.0
   r = 0.7
   RB = 0.92
   th = np.pi / 4.0
   #thh = th / 2.0
   lambda_val = 0.8
   Lz = 1
   Nc = 7
   NB = 1
   NM = 3
   compressRatio_B = 0.85
   compressRatio_M = 0.87
   Nz = 180
   c_lambda = 0.05
   RBC = RB+c_lambda

   # Calculate coordinates based on the formulas
   dx = r * np.cos(th)
   dy = r * np.sin(th)
   dxB = RB * np.cos(th)
   dyB = RB * np.sin(th)
   Dx = R * np.cos(th)
   Dy = R * np.sin(th)
   Dyx = np.sqrt((dx + lambda_val*R)**2 + dy**2) - lambda_val*R
   dxBC = RB * np.cos(th)
   dyBC = RB * np.sin(th)
   
   points_aux = np.array([
            [0, 0],
            [lambda_val * R, 0],
            [0, -lambda_val * R],
            [-lambda_val * R, 0],
            [0, lambda_val * R]
   ])
   naux = len(points_aux)

   # Block vertices
   if (half):
    points_block = np.array([
        [dx, dy],
        [dx, -dy],
        [0.0, -Dyx],
        [0.0, Dyx],
        [dxB, dyB],
        [dxB, -dyB],
        [0.0, -RB],
        [0.0, RB],
        [Dx, Dy],
        [Dx, -Dy],
        [0.0, -R],
        [0.0, R]
    ])
   else:
    points_block = np.array([
      [dx, dy],
      [dx, -dy],
      [-dx, -dy],
      [-dx, dy],
      [dxB, dyB],
      [dxB, -dyB],
      [-dxB, -dyB],
      [-dxB, dyB],
      [Dx, Dy],
      [Dx, -Dy],
      [-Dx, -Dy],
      [-Dx, Dy]
   ])
    nblk = len(points_block)

   points = np.concatenate([points_aux, points_block], axis=0)

   # Circle connections (correspond to the given circle definitions)
   if (half):
    """circles = [
        [5, 3, 6],
        [6, 4, 7],
        [7, 5, 8],
        [13, 1, 10],
        [10, 1, 11],
        [11, 1, 12],
        [17, 1, 14],
        [14, 1, 15],
        [15, 1, 16]
    ]"""
    circles = [
        [5      , 3, naux+1],
        [naux+1 , 4, naux+2],
        [naux+2 , 5, naux+3],
        [naux+8 , 6, naux+5],
        [naux+5 , 6, naux+6],
        [naux+6 , 6, naux+7],
        [naux+12, 1, naux+9],
        [naux+9 , 1, naux+10],
        [naux+10, 1, naux+11]
    ]
   else:
      """circles = [
        [9, 3, 6],
        [6, 4, 7],
        [7, 5, 8],
        [9, 2, 8],
        [13, 1, 10],
        [10, 1, 11],
        [11, 1, 12],
        [13, 1, 12],
        [17, 1, 14],
        [14, 1, 15],
        [15, 1, 16],
        [17, 1, 16]
    ]"""
      circles = [
        [4, 3, 1],
        [1, 4, 2],
        [2, 5, 3],
        [4, 2, 3],
        [8, 1, 5],
        [5, 1, 6],
        [6, 1, 7],
        [8, 1, 7],
        [12, 1, 9],
        [9, 1, 10],
        [10, 1, 11],
        [12, 1, 11]
    ]

   # Line connections (correspond to the given line definitions)
   """lines = [
      [6, 10],
      [7, 11],
      [8, 12],
      [9, 13],
      [10, 14],
      [11, 15],
      [12, 16],
      [13, 17],
   ]"""
   lines = [
      [1, 5],
      [2, 6],
      [3, 7],
      [4, 8],
      [5, 9],
      [6, 10],
      [7, 11],
      [8, 12],
   ]
   if (half):
      #lines = np.concatenate([ lines, [[8, 9]] ], axis=0)
      lines = np.concatenate([ lines, [[3, 4]] ], axis=0)

   
   # Plot the points
   fig, ax = plt.subplots(figsize=(10,8))

   # Plot auxiliary points (group 1)
   ax.scatter(points_aux[:, 0], points_aux[:, 1], color='blue', label='Auxiliary Points')

   # Plot block vertices (group 2)
   ax.scatter(points_block[:, 0], points_block[:, 1], color='red', label='Block Vertices')

   # Annotate auxiliary points
   for i, point in enumerate(points_aux):
      ax.text(point[0], point[1], f'{i+1}', color='blue', fontsize=12, ha='right', va='bottom')

   # Annotate block vertices
   for i, point in enumerate(points_block):
      ax.text(point[0], point[1], f'{points_aux.shape[0] + i+1}', color='red', fontsize=12, ha='right', va='bottom')

   # Plot circles (connecting points according to the circle definitions)
   for cidx, circle in enumerate(circles):
      p1 = points[circle[0] - 1 + naux]  # Subtract 1 for 0-indexing
      c  = points[circle[1] - 1]    #center
      p2 = points[circle[2] - 1 + naux]

      cc  = c[0] + 1j*c[1]
      p1c = p1[0] + 1j*p1[1]
      p2c = p2[0] + 1j*p2[1]
      r1  = p1c - cc
      r2  = p2c - cc
      da  = np.angle(r2) - np.angle(r1)
      if (abs(da) > np.pi):
         da = 2*np.pi - abs(da)
      alp = np.linspace(0,da,101, endpoint=True)
      circ = cc + r1*np.exp(1j*alp)
      
      # Draw the circle segments
      #ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='green', linestyle='-', linewidth=1)
      ax.plot(np.real(circ), np.imag(circ), color='green', linestyle='-', linewidth=1)
         
      # Add annotation at the midpoint
      ax.text(np.real(circ[50]), np.imag(circ[50]), f'{cidx+1}', color='green', fontsize=12, ha='left', va='bottom')

   # Plot lines (connecting points according to the line definitions)
   for lidx, line in enumerate(lines):
      p1 = points[line[0] - 1 + naux]  # Subtract 1 for 0-indexing
      p2 = points[line[1] - 1 + naux]
      
      # Draw the line connecting the two points
      ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='purple', linestyle='-', linewidth=1)

      # Calculate the midpoint for annotation
      midpoint = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2]
      
      # Add annotation at the midpoint
      ax.text(midpoint[0], midpoint[1], f'{cidx+1+lidx+1}', color='purple', fontsize=12, ha='left', va='bottom')

   # Labels and title
   ax.set_xlabel('X')
   ax.set_ylabel('Y')
   ax.set_title('Plot of Points')
   ax.legend()

   # Display the plot
   plt.axis('equal')
   plt.show()