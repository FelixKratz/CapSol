import cv2
import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt

class Mask:
    def __init__(self, x_0, x_1, y_0, y_1):
        self.x_0, self.x_1, self.y_0, self.y_1 = x_0, x_1, y_0, y_1
        if self.y_0 <= 0:
            self.y_0 = 1
        if self.y_1 <= 0:
            self.y_1 = 1
        if self.x_0 <= 0:
            self.x_0 = 1
        if self.x_1 <= 0:
            self.x_1 = 1

def plotData(ax, y, x=None, xlabel=None, ylabel=None):
    if x is None:
        ax.plot(np.array(y), marker="x")
    else:
        ax.plot(x, y)

    if xlabel is not None:
        ax.set_xlabel(xlabel)

    if ylabel is not None:
        ax.set_ylabel(ylabel)


def saveShapeToFile(file, shape):
    with open(file, "w") as f:
        for r,z in zip(shape[0].T[0], shape[0].T[1]):
            f.write(str(r) + " " + str(z) + "\n")

# Helper function to plot the resulting shapes
def plotShape(ax, shape):
    plotData(ax, x=shape[0].T[0], y=shape[0].T[1])

    ax.set_aspect("equal")
    ax.set_xlabel(r"$r / a$")
    ax.set_ylabel(r"$z / a$")

# Helper function to plot the resulting shapes
def plotShapes(ax, shapes):
    for i in range(len(shapes)):
        plotData(ax, x=shapes[i][0].T[0], y=shapes[i][0].T[1])

    ax.set_aspect("equal")
    ax.set_xlabel(r"$r / a$")
    ax.set_ylabel(r"$z / a$")

def area(shape):
    r = shape[0].T[0][1:]
    z = shape[0].T[1][1:]
    dr = r - (shape[0].T[0])[:-1]
    dz = z - (shape[0].T[1])[:-1]

    area_l = 2*np.pi * np.trapz(x=z, y=r * np.sqrt(1 + (dr / dz)**2))

    r = shape[1].T[0][1:]
    z = shape[1].T[1][1:]
    dr = r - (shape[1].T[0])[:-1]
    dz = z - (shape[1].T[1])[:-1]

    area_r = 2*np.pi * np.trapz(x=z, y=r * np.sqrt(1 + (dr / dz)**2))
    return np.array([area_l, area_r])

# Processes an array of image files into representation in parallel
def loadImages(files, flip=False, samples=250, capillary_region = 1 / 3, mask=Mask(0, 0, 0, 0)):
    return np.array(Parallel(n_jobs=multiprocessing.cpu_count())(delayed(loadImage)(file, flip, samples, capillary_region, mask) for file in files))

# Core logic to convert pixel data to the representation
def loadImage(file, flip=False, samples=250, capillary_region = 1 / 3, mask=Mask(0, 0, 0, 0), verbose=False):
    edges, grey, image_width, image_heigth = calculateEdges(file, flip=flip, mask=mask, verbose=verbose)
    cap_coords = getCapillaryViaCorners(grey, image_heigth, image_width, capillary_region)

    if (cap_coords is None):
        print("No capillary found via Corners... trying Lines")
        cap_coords = getCapillaryViaLines(edges, image_heigth, capillary_region)
    if (cap_coords is None):
        print("No capillary found via Lines...")
        return 0,0,0,0

    if verbose:
        plt.plot(cap_coords[0][0] + mask.x_0, cap_coords[0][1] + mask.y_0, marker=".", linestyle="None", color="red")
        plt.plot(cap_coords[1][0] + mask.x_0, cap_coords[1][1] + mask.y_0, marker=".", linestyle="None", color="blue")

    left_edges, right_edges, symmetry = splitAtSymmetry(edges,cap_coords)


    interp_r_s_r, interp_z_s_r,_,_,_, _, L_r, a_r = transformEdgeMatrixToInterpolation(right_edges)
    interp_r_s_l, interp_z_s_l,_, _,_,_, L_l, a_l = transformEdgeMatrixToInterpolation(left_edges)

    r_l,z_l, _ = getRepresentation(interp_r_s_l, interp_z_s_l, L_l, samples)
    r_r,z_r, _ = getRepresentation(interp_r_s_r, interp_z_s_r, L_r, samples)

    if verbose:
        plt.plot(r_r*a_r + symmetry + mask.x_0, -z_r*a_r + mask.y_0 + int(cap_coords[1][1]))
        plt.plot(-r_l*a_l + symmetry + mask.x_0, -z_l*a_l + mask.y_0 + int(cap_coords[0][1]))

    # Hacky solution to create nested numpy arrays with different sizes
    arr = np.array([[1,2,3],[4,5]], dtype=object)
    arr[0] = np.array([r_l, z_l]).T
    arr[1] = np.array([r_r, z_r]).T
    return arr

# Get the edges from pixel data
def calculateEdges(file, flip, mask, verbose):
    if isinstance(file, str):
        img = cv2.imread(file)
    else:
        img = file

    if flip:
        img = cv2.flip(img, 0)

    grey = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    grey = grey[mask.y_0 - 1:-mask.y_1]
    grey = grey.T[mask.x_0:-mask.x_1].T
    grey[:1] = np.ones(grey[:1].shape)
    edges = cv2.Canny(grey,185,200)

    if verbose:
        img[0:mask.y_0, ::] = (0,0,0)
        img[-mask.y_1:, ::] = (0,0,0)
        img[:: , 0:mask.x_0] = (0,0,0)
        img[:: , -mask.x_1:] = (0,0,0)
        plt.imshow(img)

    image_width = len(edges[0])
    image_heigth = len(edges.T[0])
    return edges, grey, image_width, image_heigth

# Capillary detection via line detection
def getCapillaryViaLines(edges, image_heigth, capillary_region):
    cap_coords = []
    lines = getLines(edges)
    for line in lines:
        x1, y1, _, y2 = line[0]

        if (y2 < image_heigth * capillary_region):
            cap_coords.append([x1, y1])

    if len(cap_coords) < 2:
        return None

    if (cap_coords[0][0] < cap_coords[1][0]):
        cap_coords.reverse()

    cap_coords = cap_coords[:2]
    cap_coords.reverse()
    return cap_coords

# Capillary detection via corner detection
def getCapillaryViaCorners(grey, image_heigth, image_width, capillary_region):
    cap_coords = []
    corners = getCorners(grey)
    for corner in corners:
        if corner[1] < image_heigth * capillary_region \
            and corner[0] > (1 / 2 - capillary_region) * image_width \
            and corner[0] < (1 / 2 + capillary_region) * image_width:
            cap_coords.append([corner[0], corner[1]])

    if len(cap_coords) < 2:
        return None

    if (cap_coords[0][0] < cap_coords[1][0]):
        cap_coords.reverse()

    cap_coords = cap_coords[:2]
    cap_coords.reverse()
    return cap_coords

def getLines(edges):
    return cv2.HoughLinesP(image=edges,rho=1,theta=1, threshold=30,
                            lines=np.array([]), minLineLength=15,
                            maxLineGap=3)
def getCorners(grey):
    dst = cv2.cornerHarris(grey,2,7,0.04)
    dst = cv2.dilate(dst,None)
    _, dst = cv2.threshold(dst,0.01*dst.max(),255,0)
    dst = np.uint8(dst)

    _, _, _, centroids = cv2.connectedComponentsWithStats(dst)
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.001)
    return cv2.cornerSubPix(grey,np.float32(centroids),(7,7),(-1,-1),criteria)

# Splits the image at symmetry
def splitAtSymmetry(edges, capillary_coordinates):
    symmetry = (capillary_coordinates[1][0] + capillary_coordinates[0][0]) / 2.

    # Proper rounding has to be performed here
    right_edges = edges.T[int(symmetry):]
    right_edges = right_edges.T
    right_edges = right_edges[int(capillary_coordinates[1][1]):]
    right_edges[0:2, ::] = 0
    right_edges[0, int(capillary_coordinates[1][0]) - int(symmetry)] = 255

    left_edges = edges.T[:int(symmetry) + 1]
    left_edges = np.flip(left_edges.T, 1)
    left_edges = left_edges[int(capillary_coordinates[0][1]):]
    left_edges[0:2, ::] = 0
    left_edges[0, int(symmetry) - int(capillary_coordinates[0][0])] = 255

    return left_edges, right_edges, symmetry

# Filter dirt or dust pixels from the shape by strictly enforcing a continous contour
def proximityFilter(raw_shape_z, raw_shape_r, proximity_threshold, symmetry_threshold):
    prox_shape_z, prox_shape_r = [], []
    found_apex = False
    for i in range(0, len(raw_shape_r)):
        if not found_apex and raw_shape_r[i] > symmetry_threshold:
            continue
        elif not found_apex:
            found_apex = True
            prox_shape_z.append(raw_shape_z[i])
            prox_shape_r.append(raw_shape_r[i])
            continue

        if (raw_shape_z[i] - prox_shape_z[-1]) * (raw_shape_z[i] - prox_shape_z[-1]) + (raw_shape_r[i] - prox_shape_r[-1]) * (raw_shape_r[i] - prox_shape_r[-1]) < proximity_threshold:
            prox_shape_z.append(raw_shape_z[i])
            prox_shape_r.append(raw_shape_r[i])
    return np.array(prox_shape_z), np.array(prox_shape_r)

def getRawShapeFromEdgeMatrix(edgeMatrix):
    raw_shape_z, raw_shape_r = np.where(np.fliplr(edgeMatrix) > 0)
    raw_shape_z, raw_shape_r = -np.flipud(raw_shape_z), len(edgeMatrix[0]) - 1 - np.flipud(raw_shape_r)

    raw_shape_z, raw_shape_r = proximityFilter(raw_shape_z, raw_shape_r, proximity_threshold=50, symmetry_threshold=15) #50 15
    return raw_shape_z, raw_shape_r

# The edge matrix is used to create a cubic interpolation of the shape
# from which *any* data representation can be gathered
def transformEdgeMatrixToInterpolation(edgeMatrix):
    # Calculate shape vectors from edge matrix
    raw_shape_z, raw_shape_r = getRawShapeFromEdgeMatrix(edgeMatrix)

    #Smoothe shape array to get subpixel precision
    savgol_shape_r = savgol_filter(raw_shape_r, 7, 2, mode="nearest") #5
    savgol_shape_z = savgol_filter(raw_shape_z, 7, 2, mode="nearest") #5

    #Scale by capillary width
    a = savgol_shape_r[-1]*2.
    savgol_shape_r = list(savgol_shape_r / a)
    savgol_shape_z = list(savgol_shape_z / a)

    z_apex = savgol_shape_z[0]

    if (savgol_shape_z[-1] != 0.):
        savgol_shape_r.append(0.5)
        savgol_shape_z.append(0.)
    else:
        savgol_shape_r[-1] = 0.5

    #Append the exact initial condition if not in list:
    if (savgol_shape_z[0] != z_apex):
        savgol_shape_r.insert(0, 0.)
        savgol_shape_z.insert(0, z_apex)
    else:
        savgol_shape_r[0] = 0.

    #Calculate deformed arc-length from smoothed shape arrays
    savgol_shape_r = np.array(savgol_shape_r)
    savgol_shape_z = np.array(savgol_shape_z)
    dx = savgol_shape_r[1::] - savgol_shape_r[:-1:]
    dy = savgol_shape_z[1::] - savgol_shape_z[:-1:]
    s = np.insert(np.cumsum(np.sqrt(dx**2 + dy**2)), 0, 0.)
    L = s[-1]

    # TODO: Filter s for duplicates before interpolating

    #Create interpolated functions r(s), z(s)
    interp_r_s = interpolate.interp1d(s, savgol_shape_r, kind='cubic')
    interp_z_s = interpolate.interp1d(s, savgol_shape_z, kind='cubic')

    return interp_r_s, interp_z_s, raw_shape_r, raw_shape_z, savgol_shape_r, savgol_shape_z, L, a

# The representation is simply calculated from the
# cubic interpolation objects and the shape count
def getRepresentation(interp_r_s, interp_z_s, L, fixedCount):
    s, r, z = [],[],[]
    if (fixedCount <= 1):
        step = 1e-2
        fixedCount = int(L / step)
        s = [step*i for i in range(fixedCount)]
        r, z = interp_r_s(s), interp_z_s(s)
    else:
        step = L / (fixedCount - 1)
        s = [step*i for i in range(fixedCount)]
        r, z = np.append(interp_r_s(s[:-1]), 0.5), np.append(interp_z_s(s[:-1]),0.)
    return r,z,s
