import numpy as np

def read_structured_vtk(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # Find DIMENSIONS line
    dim_line = next(line for line in lines if line.startswith("DIMENSIONS"))
    _, nx, ny, nz = dim_line.split()
    nx, ny, nz = int(nx), int(ny), int(nz)

    # Find POINTS section
    start_idx = next(i for i, line in enumerate(lines) if line.startswith("POINTS")) + 1

    # Read all points
    coords = []
    for line in lines[start_idx:]:
        parts = line.strip().split()
        if len(parts) == 3:
            coords.append([float(parts[0]), float(parts[1]), float(parts[2])])

    coords = np.array(coords)
    
    # Reshape to (nx, ny, nz, 3)
    coords = coords.reshape((nz, ny, nx, 3))  # note VTK order (x varies fastest)
    coords = np.transpose(coords, (2, 1, 0, 3))  # reorder to match your write loop order

    # Separate into X, Y, Z
    X = coords[..., 0]
    Y = coords[..., 1]
    Z = coords[..., 2]
    vector_list = coords.reshape(-1, 3).tolist()

        # Create vector list where y varies fastest, then x, then z
    vector_list = []
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                vector_list.append([X[i, j, k], Y[i, j, k], Z[i, j, k]])

    return vector_list, nx,ny,nz