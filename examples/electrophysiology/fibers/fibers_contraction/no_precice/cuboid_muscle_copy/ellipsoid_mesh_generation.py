import math
import matplotlib.pyplot as plt

import numpy as np

def plot_2d_projection(points):
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]

    plt.figure()
    plt.scatter(x_vals, y_vals, c='red', marker='o')
    
    # Add index labels
    for i, (x, y) in enumerate(zip(x_vals, y_vals)):
        plt.text(x, y, str(i), fontsize=9, ha='right', va='bottom', color='blue')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D Projection (Z constant)')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

def plot_3d_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Unpack x, y, z coordinates
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]
    z_vals = [p[2] for p in points]

    ax.scatter(x_vals, y_vals, z_vals, c='blue', marker='o')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D Point Visualization')
    plt.show()

# x^2/c^2 + y^2/c^2 + z^2/a^2 = 1
# 2*r^2/c^2 + z^2/a^2 = 1

# assume el_x is even, and el_x = el_y
def mesh_nodes(el_x,el_y, el_z,c,a,amin,amax):

    if el_x != el_y:
        print("Only el_x = el_y supported")

    nodes = []
    z_points = [(amax-amin)/el_z*i+amin for i in range(el_z+1)]
    print(z_points)
    for z in z_points:
        # r^2 = c²/2*(1-z^2/a^2)
        rz=np.sqrt(c**2/2*(1-z**2/a**2))
        circular_nodes = create_points(rz,el_x,z)
        for node in circular_nodes:
            nodes.append(node)


    return nodes

def compute_outer_points(el):    
    return (el + 1)*4-4

def create_points(r,el,z):
    
    nodes = []
    circular_matrix = []
    circular_matrix = circle_points(r,(el+1)*4-4,z)
    print(circular_matrix)
    
    xpoints = el +1

    if el == 2:
        # y = 0
        for i in range(xpoints):
            nodes.append(circular_matrix[i])
        # y = 1
        nodes.append(circular_matrix[-1])
        nodes.append([0,0,z])
        nodes.append(circular_matrix[i+1])
        # y = 2
        for i in range(xpoints):
            nodes.append(circular_matrix[-2-i])     
    
    if el == 4:
        xpoints_small = el-1
        circular_matrix_small = circle_points(r/2,(el-1)*4-4,z)
        
        # y = 0
        for i in range(xpoints):
            nodes.append(circular_matrix[i])
        # y = 1
        nodes.append(circular_matrix[-1])
        for j in range(xpoints_small):
            nodes.append(circular_matrix_small[j])
        nodes.append(circular_matrix[i+1])
        # y = 2
        nodes.append(circular_matrix[-2])
        nodes.append(circular_matrix_small[-1])
        nodes.append([0,0,z])
        nodes.append(circular_matrix_small[j+1])
        nodes.append(circular_matrix[i+2])

        # y = 3
        nodes.append(circular_matrix[-3])
        for j in range(xpoints_small):
            nodes.append(circular_matrix_small[-2-j])
        nodes.append(circular_matrix[i+3])
        # y = 4
        for i in range(xpoints):
            nodes.append(circular_matrix[-4-i]) 

    return nodes 

def circle_points(r, points, z):
    result = []
    for k in range(points):
        theta = math.pi + 2 * math.pi * k / points  # Start at π (180°)
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        result.append([x, y, z])
    return result



            

def apply_offset(nodes,offset):
    nodes_with_offset = []
    for node in nodes:
        nodes_with_offset.append(node + offset)
    return nodes_with_offset


################## 

#plot_2d_projection(create_points(2,4,0))