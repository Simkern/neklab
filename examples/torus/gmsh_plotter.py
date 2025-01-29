import sys, re
import numpy as np
import matplotlib.pyplot as plt
from read_2d_data import read_fields
from plot_2d_data import plot_2d_mesh

def compute_bisector(xt, yt, xb, yb):
   xm = (xt + xb)/2.0
   ym = (yt - yb)/2.0
   mr  = np.array([ xm, ym ])
   ml  = np.array([ -xm, ym ])
   # connector
   dx = xb - xt
   dy = - yb - yt
   norm = np.sqrt(dx**2 + dy**2)
   # bisector
   ux  = dy/norm  # points in neg x dir
   uy  = -dx/norm
   ul = np.array([-ux, uy])
   ur = np.array([ ux, uy])
   return ul, ur, norm

def make_circle(circle, points, npoints=101):
   p1 = points[circle[0] - 1]  # Subtract 1 for 0-indexing
   c  = points[circle[1] - 1]    #center
   p2 = points[circle[2] - 1]

   cc  = c[0] + 1j*c[1]
   p1c = p1[0] + 1j*p1[1]
   p2c = p2[0] + 1j*p2[1]
   r1  = p1c - cc
   r2  = p2c - cc
   da  = np.angle(r2) - np.angle(r1)
   if (abs(da) > np.pi):
      da = 2*np.pi - abs(da)
   alp = np.linspace(0,da,npoints, endpoint=True)
   circ = cc + r1*np.exp(1j*alp)
   #print(f'circle from {np.real(circ)[0]},{np.imag(circ)[1]} to {np.real(circ)[-1]},{np.imag(circ)[-1]}')
   return np.real(circ), np.imag(circ)

def make_line(line, points):
   p1 = points[line[0] - 1]  # Subtract 1 for 0-indexing
   p2 = points[line[1] - 1]
   #print(f'line from {p1[0]},{p1[1]} to {p2[0]},{p2[1]}')
   return [p1[0], p2[0]], [p1[1], p2[1]]

def plot_gmsh(geom, half=False, aux1=False, aux2=False, hlines=False, print_pts=False, extract_polygons=False):

   # Extract data
   R, rt, rb, RBt, RBb, tht, thb, lambda1t, lambda1b, lambda2, dyc = geom.values()

   # Calculate coordinates based on the formulas
   ra = (rt + rb)/2.0
   # aux points
   gamma = (-tht + thb)/2.0 # half angle between triple points
   lRt  = lambda1t * R
   lRb  = lambda1b * R
   lRs  = lambda2 * R
   cost = np.cos(tht)
   sint = np.sin(tht)
   cosb = np.cos(thb)
   sinb = np.sin(thb)
   cosg = np.cos(gamma)
   sing = np.sin(gamma)
   # top
   dxt  = rt * cost
   dyt  = rt * sint
   dxBt = RBt * cost
   dyBt = RBt * sint
   Dxt  = R * cost
   Dyt  = R * sint
   Dyxt = np.sqrt((dyt + lRt)**2 + dxt**2) - lRt # for the half mesh
   RBC  = np.sqrt((dyc + dyBt)**2 + dxBt**2) - dyc
   
   # bottom
   dxb  = rb * cosb
   dyb  = rb * sinb
   dxBb = RBb * cosb
   dyBb = RBb * sinb
   Dxb  = R * cosb
   Dyb  = R * sinb
   Dyxb = np.sqrt((dyb + lRb)**2 + dxb**2) - lRb # for the half mesh

   # midpoint tb1
   xm = (dxt + dxb)/2.0
   ym = (dyt - dyb)/2.0
   mr  = np.array([ xm, ym ])
   ml  = np.array([ -xm, ym ])
   # bisector
   ul, ur, norm = compute_bisector(dxt, dyt, dxb, dyb)
   h   = norm/2.0
   L   = np.sqrt((lRs + ra)**2 - h**2) # distance along bisector
   if print_pts:
      print(f'xm = {xm}')
      print(f'ym = {ym}')
      print(f'n  = {norm}')
      print(f'h  = {h}')
      print(f'ur = {ur}')
      print(f'ul = {ul}')
      print(f'L  = {L}')

   # midpoint tb2
   xmB = (dxBt + dxBb)/2.0
   ymB = (dyBt - dyBb)/2.0
   mrB  = np.array([ xmB, ymB ])
   mlB  = np.array([ -xmB, ymB ])
   # bisector
   ulB, urB, normB = compute_bisector(dxBt, dyBt, dxBb, dyBb)
   hB   = normB/2.0
   LB   = np.sqrt(RBb**2 - hB**2) # distance along bisector
   if print_pts:
      print(f'xmB = {xmB}')
      print(f'ymB = {ymB}')
      print(f'nB  = {normB}')
      print(f'hB  = {hB}')
      print(f'urB = {urB}')
      print(f'ulB = {ulB}')
      print(f'LB  = {LB}')

   points_aux2 = np.array([  mr,  ml,  mrB,  mlB ])

   points_aux = np.array([  
      [0, 0],
      ml + L*ul,
      [0, -lRt],
      mr + L*ur,
      [0, lRb],
      [0, -dyc],
      mlB + LB*ulB,
      mrB + LB*urB
   ])
   naux = len(points_aux)

   # Block vertices
   points = np.array([
      [  0.0, Dyxt], 
      [  dxt,  dyt], 
      [  dxb, -dyb], 
      [  0.0,-Dyxb], 
      [  0.0,  RBC], 
      [ dxBt, dyBt], 
      [ dxBb,-dyBb], 
      [  0.0, -RBb], 
      [  0.0,    R], 
      [  Dxt,  Dyt], 
      [  Dxb, -Dyb], 
      [  0.0,   -R]
   ])
   npts = len(points)

   points_sym = np.array([
      [ -dxb, -dyb], 
      [ -dxt,  dyt], 
      [-dxBb,-dyBb], 
      [-dxBt, dyBt], 
      [ -Dxb, -Dyb], 
      [ -Dxt,  Dyt]
   ])

   pts = np.concatenate([points, points_aux, points_sym], axis=0)

   b0 = [ 
      points_aux2[0],
      points_aux2[2]
   ]
   l0 = [ 
      points[2-1],
      points[3-1],
      points[6-1],
      points[7-1],
   ]
   
   if not half:
      b0 += [
         points_aux2[1],
         points_aux2[3]
      ]
      l0 += [
         points_sym[21-naux-npts-1],
         points_sym[22-naux-npts-1],
         points_sym[25-naux-npts-1],
         points_sym[26-naux-npts-1]
      ]

   b1 = [ points_aux[3], points_aux[6] ]
   l1 = 2*[ points_aux[3] ] + 2*[ points_aux[6] ]
   if not half:
      b1 += [ points_aux[1], points_aux[7] ]
      l1 += 2*[ points_aux[1] ] + 2*[ points_aux[7] ]

   bline0 = np.array(b0)
   hline0 = np.array(l0)
   bline1 = np.array(b1)
   hline1 = np.array(l1)

   if print_pts:
      p = -points[1]
      pc = (p[0] +1j*p[1])/np.linalg.norm(p)
      p = points[1] - points_aux[3]
      p1 = (p[0] +1j*p[1])/np.linalg.norm(p)
      p = points[1] - points_aux[2]
      p2 = (p[0] +1j*p[1])/np.linalg.norm(p)
      print(f'pc = {pc}')
      print(f'p1 = {p1}')
      print(f'p2 = {p2}')
      q1 = np.real(np.dot(pc,p1))
      q2 = np.real(np.dot(p2,pc))
      print(f'cosa = {q1}')
      print(f'cosb = {q2}')
      print(f'a = {np.arccos(q1)*360/(2*np.pi)}')
      print(f'b = {180-np.arccos(q2)*360/(2*np.pi)}')

      p = -points[2]
      pc = (p[0] +1j*p[1])/np.linalg.norm(p)
      p = points[2] - points_aux[3]
      p1 = (p[0] +1j*p[1])/np.linalg.norm(p)
      p = points[2] - points_aux[4]
      p2 = (p[0] +1j*p[1])/np.linalg.norm(p)
      print(f'pc = {pc}')
      print(f'p1 = {p1}')
      print(f'p2 = {p2}')
      q1 = np.real(np.dot(pc,p1))
      q2 = np.real(np.dot(p2,pc))
      print(f'cosa = {q1}')
      print(f'cosb = {q2}')
      print(f'a = {np.arccos(q1)*360/(2*np.pi)}')
      print(f'b = {180-np.arccos(q2)*360/(2*np.pi)}')

   # Circle connections (correspond to the given circle definitions)
   circles = [
      # inner ring
      [1, 15, 2],
      [2, 16, 3],
      [3, 17, 4],
      # middle ring
      [5, 18, 6],
      [6, 20, 7],
      [7, 13, 8],
      # outer ring
      [ 9, 13, 10],
      [10, 13, 11],
      [11, 13, 12]
   ]
   if not half:
      circles += [
      # inner ring
         [ 4, 17, 21],
         [22, 14, 21], # flipped order for plot
         [22, 15,  1],
         [ 8, 13, 23],
         [24, 19, 23], # flipped order for plot
         [24, 18,  5],
         [12, 13, 25],
         [26, 13, 25], # flipped order for plot
         [26, 13,  9]
      ]

   # Line connections (correspond to the given line definitions)
   lines = [
      [1, 5],
      [2, 6],
      [3, 7],
      [4, 8],
      [5, 9],
      [6, 10],
      [7, 11],
      [8, 12],
      [ 4, 1]
   ]
   if not half:
      lines += [
         [ 21, 23 ],
         [ 22, 24 ],
         [ 23, 25 ],
         [ 24, 26 ]
      ]
   
   # Plot the points
   fig, ax = plt.subplots(figsize=(10,8))

   if aux1:
      # Plot auxiliary points (group 1)
      ax.scatter(points_aux[:, 0], points_aux[:, 1], color='blue', label='Auxiliary Points')
      # Annotate auxiliary points
      for i, point in enumerate(points_aux):
         ax.text(point[0], point[1], f'{i+1}', color='blue', fontsize=12, ha='right', va='bottom')
   
   # Plot block vertices (group 2)
   ax.scatter(points[:, 0], points[:, 1], color='red', label='Block Vertices')  

   # Annotate block vertices
   for i, point in enumerate(points):
      ax.text(point[0], point[1], f'{i+1}', color='red', fontsize=12, ha='right', va='bottom')

   if not half:
      # Plot block vertices (group 2)
      ax.scatter(points_sym[:, 0], points_sym[:, 1], color='red')  

      # Annotate block vertices
      for i, point in enumerate(points_sym):
         ax.text(point[0], point[1], f'{naux + points.shape[0] + i+1}', color='red', fontsize=12, ha='right', va='bottom')

   if aux2:
      # Plot auxiliary points 2
      ax.scatter(points_aux2[:, 0], points_aux2[:, 1], color='black', label='Aux2')
      for ia, point in enumerate(points_aux2):
         ax.text(point[0], point[1], f'{ia+1}', color='black', fontsize=12, ha='left', va='bottom')
   
   if hlines:
      for i, (px, py) in enumerate(zip(bline0,bline1)):
         ax.plot([px[0], py[0]], [px[1], py[1]], c = 'black', linestyle='--')
      for i, (px, py) in enumerate(zip(hline0,hline1)):
         ax.plot([px[0], py[0]], [px[1], py[1]], c = 'black', alpha=0.25)
   
   # Plot circles (connecting points according to the circle definitions)
   for cidx, circle in enumerate(circles):
      xcirc, ycirc = make_circle(circle, pts)
      
      # Draw the circle segments
      ax.plot(xcirc, ycirc, color='green', linestyle='-', linewidth=1)
         
      # Add annotation at the midpoint
      ax.text(xcirc[50], ycirc[50], f'{cidx+1}', color='green', fontsize=12, ha='left', va='bottom')

   # Plot lines (connecting points according to the line definitions)
   for lidx, line in enumerate(lines):
      xp, yp = make_line(line, pts)
      
      # Draw the line connecting the two points
      ax.plot(xp, yp, color='purple', linestyle='-', linewidth=1)

      # Calculate the midpoint for annotation
      midpoint = [np.sum(xp)/2.0, np.sum(yp)/2.0]
      
      # Add annotation at the midpoint
      ax.text(midpoint[0], midpoint[1], f'{cidx+1+lidx+1}', color='purple', fontsize=12, ha='left', va='bottom')
 
   # Labels and title
   ax.set_xlabel('X')
   ax.set_ylabel('Y')
   ax.set_title('Plot of Points')
   ax.legend()

   # Display the plot
   plt.axis('equal')