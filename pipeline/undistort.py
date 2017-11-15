import cv2
import numpy as np
import argparse
import os
def undistort(img_path, d_lambda, is_full):
    img = cv2.imread(img_path)
    rows, cols = img.shape[:2]

    r_img = np.sqrt((rows / 2.0) ** 2 + (cols / 2.0) ** 2)
    m_ind = np.indices((rows, cols), np.double)
    if (is_full <= 0):
        d = min(rows, cols)
    else:
        d = max(rows, cols)

    alpha = 4/(4+d_lambda*(d/r_img)**2)

    x_ind = alpha*(m_ind[1, :, :] - cols / 2.0) / (r_img)
    y_ind = alpha*(m_ind[0, :, :] - rows / 2.0) / (r_img)
    m_r_u = x_ind ** 2 + y_ind ** 2 + 1e-9

    m_r_d = np.divide(1 - np.sqrt(1 - 4 * d_lambda * m_r_u), 2 * d_lambda * np.sqrt(m_r_u))
    map_x = (r_img* x_ind * (1 + d_lambda * m_r_d ** 2) + cols / 2.0).astype(np.float32)
    map_y = (r_img* y_ind * (1 + d_lambda * m_r_d ** 2) + rows / 2.0).astype(np.float32)

    dst = cv2.remap(img, map_x, map_y, cv2.INTER_LINEAR)
    name_with_extension = os.path.basename(img_path)
    input_directory = os.path.split(os.path.dirname(img_path))[1]
    name = os.path.splitext(name_with_extension)[0]
    if (is_full <= 0):
        out_directory ='./undistorted_images/{}/'.format(input_directory)
    else:
        out_directory ='./full_undistorted_images/{}/'.format(input_directory)

    if (not os.path.exists(out_directory )):
        os.mkdir(out_directory)

    undist_img_path = '{}{}_undistorted_k0_is_{:0.6f}.jpg'.format(out_directory, name, d_lambda)
    print('Output file path: ' + undist_img_path)
    cv2.imwrite(undist_img_path, dst)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Undistort image')
    parser.add_argument('--coeff', dest='distort_coeff', type=str, required=True)
    parser.add_argument('--img', dest='img_path', type=str, required=True)
    parser.add_argument('--full', dest='is_full_out', type=int, default=0)
    args = parser.parse_args()
    print('Undistortion with parameter: ' + str(args.distort_coeff))
    undistort(args.img_path, float(args.distort_coeff), args.is_full_out)

