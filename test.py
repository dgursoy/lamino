import lamino
import numpy as np
import dxchange
import tomopy

# Define object's values
obj = dxchange.read_tiff('/Users/dgursoy/Data/Chip/lamino/beta.tiff')[1:71, 1:71, 1:71]
obj = np.pad(obj, ((15, 15), (15, 15), (15, 15)), mode='constant')
dxchange.write_tiff(obj)

# Define object's grid coordinates
gridx = np.linspace(-0.5, 0.5, obj.shape[0]+1, dtype='float32')
gridy = np.linspace(-0.5, 0.5, obj.shape[1]+1, dtype='float32')
gridz = np.linspace(-0.5, 0.5, obj.shape[2]+1, dtype='float32')

# Define tilt and rotation angles
NUMPRJ = 100
phi = np.linspace(45, 45, NUMPRJ, dtype='float32') * np.pi / 180. # Tilt
theta = np.linspace(0, 360, NUMPRJ, dtype='float32') * np.pi / 180. # Rotation

# Define detector pixellation
detgridx = np.linspace(-0.5, 0.5, 100, dtype='float32')
detgridy = np.linspace(-0.5, 0.5, 100, dtype='float32')

# Calculate projections
prj = lamino.project(obj, gridx, gridy, gridz, phi, theta, detgridx, detgridy)
dxchange.write_tiff(prj)

# Reconstruction grid
rgridx = np.linspace(-0.5, 0.5, 101, dtype='float32')
rgridy = np.linspace(-0.5, 0.5, 101, dtype='float32')
rgridz = np.linspace(-0.5, 0.5, 101, dtype='float32')

# Calculate reconstruction
rec = lamino.mlem(prj, phi, theta, detgridx, detgridy, rgridx, rgridy, rgridz, iters=10)
dxchange.write_tiff(rec)
