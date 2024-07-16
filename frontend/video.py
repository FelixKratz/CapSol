import cv2
import numpy as np
from image import loadImage, Mask

def loadVideo(path, flip=False, capillary_region= 1 / 3, samples=250, mask=Mask(0, 0, 0 ,0), firstFrame=0, finalFrame=-1, verbose=None):
    vid = cv2.VideoCapture(path)
    
    result = []
    vb = False
    success = True
    counter = 0
    while(vid.isOpened()):
        counter += 1
        success, frame = vid.read()
        if not success:
            break

        if counter < firstFrame - 1:
            continue

        if verbose is not None and counter == verbose:
            vb = True
        elif verbose is not None:
            continue

        r = loadImage(frame,
                      flip=flip,
                      samples=samples,
                      capillary_region=capillary_region,
                      mask=mask,
                      verbose=vb)

        if vb:
            return r

        result.append(r)
        if finalFrame > 0 and finalFrame == counter - 1:
            break

    vid.release()
    return np.array(result)
