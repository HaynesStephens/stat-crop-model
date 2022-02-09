import numpy as np


def fillIn(K):
    if K.shape[0] == 34:
        print('CHILL.')

    elif K.shape[0] == 10:
        fa      = np.zeros((34, K.shape[1], K.shape[2]))
        fa_is   = [0, 1, 2, 5, 6, 9, 15, 16, 19, 25]
        for i in range(10):
            fa_i     = fa_is[i]
            fa[fa_i, :, :] = K[i, :, :]

    elif K.shape[0] == 19:
        fa      = np.zeros((34, K.shape[1], K.shape[2]))
        fa_is   = [0, 1, 2, 4, 5, 6, 8, 9, 11, 14, 15, 16, 18, 19, 21, 24, 25, 27, 30]
        for i in range(19):
            fa_i     = fa_is[i]
            fa[fa_i, :, :] = K[i, :, :]

    elif K.shape[0] == 20:
        fa = np.zeros((34, K.shape[1], K.shape[2]))
        fa_is = [0, 1, 2, 3, 5, 6, 7, 9, 10, 12, 15, 16, 17, 19, 20, 22, 25, 26, 28, 31]
        for i in range(20):
            fa_i = fa_is[i]
            fa[fa_i, :, :] = K[i, :, :]

    else:
        print('FUCK!')
    return fa


def transformAnom(K):
    c = 360
    w = 1
    n = 200
    k = {}
    out_array = np.zeros(K.shape)
    if K.shape[0] == 34:
        # get k values in a dictionary for ease
        for i in range(K.shape[0]):
            k[i] = K[i, :, :]

        # reorder parameters now using a anomaly framework
        out_array[0, :, :] = (
                        k[0] + (k[1] * c) + (k[5] * c ** 2) + (k[15] * c ** 3) + (k[4] * n) + (k[8] * c * n) +
                        (k[18] * n * c ** 2) + (k[14] * n ** 2) + (k[24] * c * n ** 2) + (k[3] * w) + (k[7] * w * c) +
                        (k[17] * w * c ** 2) + (k[13] * n * w) + (k[23] * c * n * w) + (k[33] * w * n ** 2) +
                        (k[12] * w ** 2) + (k[22] * c * w ** 2) + (k[32] * n * w ** 2) + (k[31] * w ** 3)
        )
        out_array[1, :, :] = (
                        k[1] + (k[5] * 2 * c) + (k[15] * 3 * c ** 2) + (k[8] * n) + (k[18] * 2 * c * n) +
                        (k[24] * n ** 2) + (k[7] * w) + (k[17] * 2 * c * w) + (k[23] * n * w) + (k[22] * w ** 2)
        )
        out_array[2, :, :] = (
                        k[2] + (k[6] * c) + (k[16] * c ** 2) + (k[11] * n) + (k[21] * c * n) + (k[30] * n ** 2) +
                        (k[10] * w) + (k[20] * c * w) + (k[29] * n * w) + (k[28] * w ** 2)
        )
        out_array[3, :, :] = (
                        k[3] + (k[7] * c) + (k[17] * c ** 2) + (k[13] * n) + (k[23] * c * n) + (k[33] * n ** 2) +
                        (k[12] * 2 * w) + (k[22] * 2 * c * w) + (k[32] * 2 * n * w) + (k[31] * 3 * w ** 2)
        )
        out_array[4, :, :] = (
                        k[4] + (k[8] * c) + (k[18] * c ** 2) + (k[14] * 2 * n) + (k[24] * 2 * c * n) + (k[13] * w) +
                        (k[23] * c * w) + (k[33] * 2 * n * w) + (k[32] * w ** 2)
        )
        out_array[5, :, :] = (
                        k[5] + (k[15] * 3 * c) + (k[18] * n) + (k[17] * w)
        )
        out_array[6, :, :] = (
                        k[6] + (k[16] * 2 * c) + (k[21] * n) + (k[20] * w)
        )
        out_array[7, :, :] = (
                        k[7] + (k[17] * 2 * c) + (k[23] * n) + (k[22] * 2 * w)
        )
        out_array[8, :, :] = (
                        k[8] + (k[18] * 2 * c) + (k[24] * 2 * n) + (k[23] * w)
        )
        out_array[9, :, :] = (
                        k[9] + (k[19] * c) + (k[27] * n) + (k[26] * w)
        )
        out_array[10, :, :] = (
                        k[10] + (k[20] * c) + (k[29] * n) + (k[28] * 2 * w)
        )
        out_array[11, :, :] = (
                        k[11] + (k[21] * c) + (k[30] * 2 * n) + (k[29] * w)
        )
        out_array[12, :, :] = (
                        k[12] + (k[22] * c) + (k[32] * n) + (k[31] * 3 * w)
        )
        out_array[13, :, :] = (
                        k[13] + (k[23] * c) + (k[33] * 2 * n) + (k[32] * 2 * w)
        )
        out_array[14, :, :] = (
                        k[14] + (k[24] * c) + (k[33] * w)
        )
        # no change in parameters beyond the 14th index, so set the rest equal to K
        out_array[15:, :, :] = K[15:, :, :]

    else: # if the shape of the array isn't 34, then expand and fill it with zeros, then transform it using anomalies
        print('N:', K.shape[0])
        fa = fillIn(K) # fill in K to make a 34 parameter array
        out_array = transformAnom(fa) # rerun transformAnom on the new, 34 parameter array

    return out_array


def runAnoms():
    import glob
    filedir = '/project2/moyer/Haynes/emulator_params/'
    savedir = '/project2/moyer/Haynes/emulator_params/anoms34/'
    file_list = glob.glob(filedir + 'LPJmL_soy*.npy')
    for i in file_list:
        name = i.split('/')[-1].split('.')[0] # get filename
        print(name)
        K = np.load(i) # load parameters
        out_array = transformAnom(K)
        if out_array is None:
            print(i)
        else:
            np.save(savedir + name + '.npy', out_array)


def emulate(K, C, T, W, N):
    if K.shape[0] == 34:
        Y = (K[0, :, :] + K[1, :, :] * C + K[2, :, :] * T + K[3, :, :] * W + K[4, :, :] * N +
             K[5, :, :] * C ** 2 + K[6, :, :] * C * T + K[7, :, :] * C * W + K[8, :, :] * C * N +
             K[9, :, :] * T ** 2 + K[10, :, :] * T * W + K[11, :, :] * T * N +
             K[12, :, :] * W ** 2 + K[13, :, :] * W * N +
             K[14, :, :] * N ** 2 +
             K[15, :, :] * C ** 3 + K[16, :, :] * C ** 2 * T + K[17, :, :] * C ** 2 * W + K[18, :, :] * C ** 2 * N +
             K[19, :, :] * C * T ** 2 + K[20, :, :] * C * T * W + K[21, :, :] * C * T * N +
             K[22, :, :] * C * W ** 2 + K[23, :, :] * C * W * N +
             K[24, :, :] * C * N ** 2 +
             K[25, :, :] * T ** 3 + K[26, :, :] * T ** 2 * W + K[27, :, :] * T ** 2 * N +
             K[28, :, :] * T * W ** 2 + K[29, :, :] * T * W * N + K[30, :, :] * T * N ** 2 +
             K[31, :, :] * W ** 3 + K[32, :, :] * W ** 2 * N +
             K[33, :, :] * W * N ** 2)
             # No N**3 term??

    elif K.shape[0] == 20:
        Y = (K[0, :, :] + K[1, :, :] * C + K[2, :, :] * T + K[3, :, :] * W +
             K[4, :, :] * C ** 2 + K[5, :, :] * C * T + K[6, :, :] * C * W +
             K[7, :, :] * T ** 2 + K[8, :, :] * T * W +
             K[9, :, :] * W ** 2 +
             K[10, :, :] * C ** 3 + K[11, :, :] * C ** 2 * T + K[12, :, :] * C ** 2 * W +
             K[13, :, :] * C * T ** 2 + K[14, :, :] * C * T * W +
             K[15, :, :] * C * W ** 2 +
             K[16, :, :] * T ** 3 + K[17, :, :] * T ** 2 * W + K[18, :, :] * T * W ** 2 +
             K[19, :, :] * W ** 3)

    elif K.shape[0] == 19:
        Y = (K[0, :, :] + K[1, :, :] * C + K[2, :, :] * T + K[3, :, :] * N +
             K[4, :, :] * C ** 2 + K[5, :, :] * C * T + K[6, :, :] * C * N +
             K[7, :, :] * T ** 2 + K[8, :, :] * T * N +
             K[9, :, :] * N ** 2 +
             K[10, :, :] * C ** 3 + K[11, :, :] * C ** 2 * T + K[12, :, :] * C ** 2 * N +
             K[13, :, :] * C * T ** 2 + K[14, :, :] * C * T * N +
             K[15, :, :] * C * N ** 2 +
             K[16, :, :] * T ** 3 + K[17, :, :] * T ** 2 * N +
             K[18, :, :] * T * N ** 2)

    elif K.shape[0] == 10:
        Y = (K[0, :, :] + K[1, :, :] * C + K[2, :, :] * T +
             K[3, :, :] * C ** 2 + K[4, :, :] * C * T +
             K[5, :, :] * T ** 2 +
             K[6, :, :] * C ** 3 + K[7, :, :] * C ** 2 * T +
             K[8, :, :] * C * T ** 2 +
             K[9, :, :] * T ** 3)

    Y = np.nan_to_num(Y)
    Y[Y < 0.01] = 0
    return (Y)
