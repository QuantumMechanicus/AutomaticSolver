#!/usr/bin/python
import cv2
import numpy as np
import argparse
import os
import pandas as pd
def undistort(img_path, f_x, f_y):
    img = cv2.imread(img_path)
    rows, cols = img.shape[:2]


    df1 = pd.read_table(f_x, delimiter=',', header=None, dtype=np.float32)
    df2 = pd.read_table(f_y, delimiter=',',header=None,  dtype=np.float32)


    map_x = df1.as_matrix()
    map_y = df2.as_matrix()
    print(map_y.shape)
    dst = cv2.remap(img, map_x, map_y, cv2.INTER_LINEAR)
    name_with_extension = os.path.basename(img_path)
    input_directory = os.path.split(os.path.dirname(img_path))[1]
    name = os.path.splitext(name_with_extension)[0]




    undist_img_path = '{}_undistorted.jpg'.format(name)
    print('Output file path: ' + undist_img_path)
    cv2.imwrite(undist_img_path, dst)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Undistort image')
    parser.add_argument('--map_x', dest='m_x', type=str, required=True)
    parser.add_argument('--map_y', dest='m_y', type=str, required=True)
    parser.add_argument('--img', dest='img_path', type=str, required=True)

    args = parser.parse_args()
    print('Undistortion')
    undistort(args.img_path, args.m_x, args.m_y)

