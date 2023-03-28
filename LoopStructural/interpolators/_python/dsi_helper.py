import numpy as np


def cg(EG: np.ndarray, neighbours: np.ndarray, elements: np.ndarray, nodes, region):
    if EG.shape[0] != neighbours.shape[0] and EG.shape[0] != elements.shape[0]:
        raise ValueError("EG and neighbours must have the same number of elements")
    if EG.shape[2] != 4:
        raise ValueError("EG must have 4 columns")
    Nc = 5  # numer of constraints shared nodes + independent
    Na = 4  # number of nodes
    Ns = Na - 1
    ne = len(neighbours)
    ncons = 0
    flag = np.zeros(ne, dtype=np.int32)
    c = np.zeros((len(neighbours) * 4, Nc))
    areas = np.zeros((len(neighbours) * 4))
    idc = np.zeros((ne * 4, 5), dtype=np.int64)
    common = np.zeros((3), dtype=np.int64)
    norm = np.zeros((3))
    shared_pts = np.zeros((3, 3))
    v1 = np.zeros(3)
    v2 = np.zeros(3)
    area = 0
    length = 0
    idl = np.zeros(4, dtype=np.int64)
    idr = np.zeros(4, dtype=np.int64)
    for e in range(ne):
        idl = elements[e, :]
        e1 = EG[e, :, :]
        flag[e] = 1
        # if not in region then skip this tetra
        if (
            region[idl[0]] == 0
            or region[idl[1]] == 0
            or region[idl[2]] == 0
            or region[idl[3]] == 0
        ):
            continue
        for n in range(4):
            neigh = neighbours[e, n]
            idr = elements[neigh, :]
            if neigh < 0:
                continue
            if flag[neigh] == 1:
                continue

            # if not in region then skip this tetra
            if (
                region[idr[0]] == 0
                or region[idr[1]] == 0
                or region[idr[2]] == 0
                or region[idr[3]] == 0
            ):
                continue

            e2 = EG[neigh, :, :]

            for i in range(Nc):
                idc[ncons, i] = -1

            i = 0
            for itr_right in range(Na):
                for itr_left in range(Na):
                    if idl[itr_left] == idr[itr_right]:
                        common[i] = idl[itr_left]
                        i += 1
            for j in range(3):
                for k in range(3):
                    shared_pts[j][k] = nodes[common[j]][k]  # common
            for i in range(3):
                v1[i] = shared_pts[0, i] - shared_pts[1, i]
                v2[i] = shared_pts[2, i] - shared_pts[1, i]
            norm[0] = v2[2] * v1[1] - v1[2] * v2[1]
            norm[1] = v1[2] * v2[0] - v1[0] * v2[2]
            norm[2] = v1[0] * v2[1] - v1[1] * v2[0]

            length = np.linalg.norm(norm)
            # we want to weight the cg by the area of the shared face
            # area of triangle is half area of parallelogram
            # https://math.stackexchange.com/questions/128991/how-to-calculate-the-area-of-a-3d-triangle
            area = (
                0.5 * length
            )  # sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2])#np.linalg.norm(norm)
            for i in range(3):
                norm[i] /= length
            for itr_left in range(Na):
                idc[ncons, itr_left] = idl[itr_left]
                for i in range(3):
                    c[ncons, itr_left] += norm[i] * e1[i][itr_left]
            next_available_position = Na
            for itr_right in range(Na):
                common_index = -1
                for itr_left in range(Na):
                    if idc[ncons, itr_left] == idr[itr_right]:
                        common_index = itr_left

                position_to_write = 0
                if common_index != -1:
                    position_to_write = common_index
                else:
                    position_to_write = 4  # next_available_position
                    next_available_position += 1
                idc[ncons, position_to_write] = idr[itr_right]
                for i in range(3):
                    c[ncons, position_to_write] -= norm[i] * e2[i][itr_right]
            areas[ncons] = area
            ncons += 1
    return idc, c, ncons, areas


def fold_cg(EG, X, neighbours, elements, nodes):
    Nc = 5  # numer of constraints shared nodes + independent
    Na = 4  # number of nodes
    Ns = Na - 1
    ne = len(neighbours)
    ncons = 0
    flag = np.zeros(ne, dtype=np.int32)
    c = np.zeros((len(neighbours) * 4, Nc))
    idc = np.zeros((ne * 4, 5), dtype=np.int64)
    common = np.zeros((3), dtype=np.int64)
    norm = np.zeros((3))
    shared_pts = np.zeros((3, 3))
    v1 = np.zeros(3)
    v2 = np.zeros(3)

    idl = np.zeros(4, dtype=np.int64)
    idr = np.zeros(4, dtype=np.int64)
    for e in range(ne):
        idl = elements[e, :]
        e1 = EG[e, :, :]
        flag[e] = 1
        Xl = X[e, :]
        for n in range(4):
            neigh = neighbours[e, n]
            idr = elements[neigh, :]
            if neigh == -1:
                continue
            if flag[neigh] == 1:
                continue
            e2 = EG[neigh, :, :]
            Xr = X[neigh, :]

            for i in range(Nc):
                idc[ncons, i] = -1
            i = 0
            for itr_right in range(Na):
                for itr_left in range(Na):
                    if idl[itr_left] == idr[itr_right]:
                        common[i] = idl[itr_left]
                        i += 1
            for j in range(3):
                for k in range(3):

                    shared_pts[j][k] = nodes[common[j]][k]  # common
            for i in range(3):
                v1[i] = shared_pts[0, i] - shared_pts[1, i]
                v2[i] = shared_pts[2, i] - shared_pts[1, i]
            norm[0] = v2[2] * v1[1] - v1[2] * v2[1]
            norm[1] = v1[2] * v2[0] - v1[0] * v2[2]
            norm[2] = v1[0] * v2[1] - v1[1] * v2[0]
            area = 0.5 * np.sqrt(
                norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]
            )  # np.linalg.norm(norm)

            i = 0
            for itr_left in range(Na):
                idc[ncons, itr_left] = idl[itr_left]
                for i in range(3):
                    c[ncons, itr_left] += Xl[i] * e1[i][itr_left] * area
            next_available_position = Na
            for itr_right in range(Na):
                common_index = -1
                for itr_left in range(Na):
                    if idc[ncons, itr_left] == idr[itr_right]:
                        common_index = itr_left
                position_to_write = 0
                if common_index != -1:
                    position_to_write = common_index
                else:
                    position_to_write = 4  # next_available_position
                    next_available_position += 1

                idc[ncons, position_to_write] = idr[itr_right]
                for i in range(3):
                    c[ncons, position_to_write] -= Xr[i] * e2[i][itr_right] * area
            ncons += 1
    return idc, c, ncons
