import os
import sys
import numpy as np
import PIL.Image as pimg

class CropJob(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def load_coordinates_file(coords_file, images_dir, cropped_dir):
    jobs = []
    with open(coords_file) as cf:
        for line in cf:
            fields = line.split()

            if len(fields) == 5:
                fields += [0]*4

            (image,
            lex, ley,
            rex, rey,
            nx, ny,
            mx, my) = [fields[0]] + map(float, fields[1:])

            jobs.append(
                CropJob(
                    image = os.path.join(images_dir, image + '.jpg'),
                    output = os.path.join(cropped_dir, image + '.jpg'),
                    left_eye = np.array((lex, ley), dtype = 'int'),
                    right_eye = np.array((rex, rey), dtype = 'int'),
                    nose = np.array((nx, ny), dtype = 'int'),
                    mouth = np.array((mx, my), dtype = 'int')
                )
            )

    return jobs

def process_crop_job(job, new_size):
    '''
    Crop an image to retain as much of the face as possible, and as little of
    the background as possible.

    Uses the proportions found at https://tinyurl.com/mfe9ams
    '''

    image = pimg.open(job.image)

    # the distance between the eyes is the width of one eye
    one_eye_width = np.linalg.norm(job.right_eye - job.left_eye) / 2

    # the face is five eyes wide
    face_width = 5 * one_eye_width

    mid_eyes = (job.left_eye + job.right_eye) / 2
    bleft = mid_eyes[0] - face_width / 2
    bright = mid_eyes[0] + face_width / 2

    mid_eyes_to_nose = np.linalg.norm(mid_eyes - job.nose)

    # the eyes are half way between the top of the head and chin, and
    # the tip of the nose is half way between the eyes and chin
    # FIXME this should really be 4 * mid_eyes_to_nose, but 8 works
    # a lot better for all images!
    face_height = 8 * mid_eyes_to_nose

    btop = mid_eyes[1] - face_height / 2
    bbottom = mid_eyes[1] + face_height / 2

    box = map(int, (bleft, btop, bright, bbottom))
    image.crop(box).resize(new_size, pimg.ANTIALIAS).save(job.output)

if __name__ == '__main__':
    images_dir = sys.argv[1]
    coords_file = sys.argv[2]
    cropped_dir = sys.argv[3]

    size = (100, 100)
    if len(sys.argv) == 5:
        size = map(int, sys.argv[4].split('x'))

    jobs = load_coordinates_file(coords_file, images_dir, cropped_dir)
    for i, j in enumerate(jobs):
        print '\rProcessing job {} / {}'.format(i + 1, len(jobs)),
        process_crop_job(j, size)
