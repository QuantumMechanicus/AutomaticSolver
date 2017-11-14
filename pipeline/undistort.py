import cv2
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Undistort image')
    parser.add_argument('--coeff', dest='distort_coeff', type=str, required=True)
    parser.add_argument('--img', dest='img_path', type=str, required=True)
    parser.add_argument('--output', dest='img_path_out', type=str, required=True)
    args = parser.parse_args()
    print(args.distort_coeff)
    img = cv2.imread(args.img_path)
    ind = 0
    d_lambda = float(args.distort_coeff)
    rows, cols = img.shape[:2]
    map_x = np.zeros((rows, cols), np.double)
    map_y = np.zeros((rows, cols), np.double)

    r_img = np.sqrt((rows / 2.0) ** 2 + (cols / 2.0) ** 2)
    m_ind = np.indices((rows, cols), np.double)
    x_ind = (m_ind[1, :, :] - cols / 2.0) / r_img
    y_ind = (m_ind[0, :, :] - rows / 2.0) / r_img
    m_r_u = x_ind ** 2 + y_ind ** 2 + 1e-9

    m_r_d = np.divide(1 - np.sqrt(1 - 4 * d_lambda * m_r_u), 2 * d_lambda * np.sqrt(m_r_u))
    map_x = (r_img * x_ind * (1 + d_lambda * m_r_d ** 2) + cols / 2.0).astype(np.float32)
    map_y = (r_img * y_ind * (1 + d_lambda * m_r_d ** 2) + rows / 2.0).astype(np.float32)

    dst = cv2.remap(img, map_x, map_y, cv2.INTER_LINEAR)
    print('Undistortion\'ve done')
    str = args.img_path_out + '_undistorted_l_is_' + str(d_lambda) + '.jpg'
    cv2.imwrite(str, dst)
