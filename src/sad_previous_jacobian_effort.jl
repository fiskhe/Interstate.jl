################################
################################
######### SAD NESSS!!! #########
################################
################################

function h_jacobian(bbox_i, points, camera, full_state)
    j = zeros(4,8)
    (left, top, right, bottom) = bbox_i
    R = camera.R
    count = 1
    # everything * f/camera.sx or sy
    for i in [left, right]
        row = zeros(8)
        if i != 0
            (pt, lwh_c) = points[i]
            # x -> 1/3
            n = pt[1]
            d = pt[3]

            config = ones(8)
            config[4:6] = lwh_c
            s = copy(full_state) .* config # shows which corner

            row[1] = dh_dx(1, n, d, R, s) # x
            row[2] = dh_dy(1, n, d, R, s) # y
            row[3] = dh_dθ(1, n, d, R, s)
            row[4] = dh_dl(1, n, d, R, s)
            row[5] = dh_dw(1, n, d, R, s)
            row[6] = dh_dh(1, n, d, R, s)
            row .* camera.focal_len/camera.sx
        end
        j[count, :] = row
        count += 1
    end
    for i in [top, bottom]
        row = zeros(8)
        if i != 0
            (pt, lwh_c) = points[i]
            # x -> 2/3
            n = pt[2]
            d = pt[3]

            config = ones(8)
            config[4:6] = lwh_c
            s = copy(full_state) .* config

            row[1] = dh_dx(2, n, d, R, s) # x
            row[2] = dh_dy(2, n, d, R, s) # y
            row[3] = dh_dθ(2, n, d, R, s)
            row[4] = dh_dl(2, n, d, R, s)
            row[5] = dh_dw(2, n, d, R, s)
            row[6] = dh_dh(2, n, d, R, s)
            row .* camera.focal_len/camera.sy
        end
        j[count, :] = row
        count += 1
    end

    # so that it is left top right bottom instead of left right top bottom
    row2 = copy(j[2, :])
    j[2, :] = j[3, :]
    j[3, :] = row2
    return j
end

function dh_dx(o, n, d, R, x)
    # to account for min/max
    if x[1] == 1.0 || x[1] == -1.0
        0
    else
        u = 3
        (d*R[o, 1] - n*R[u, 1]) / d^2
    end
end

function dh_dy(o, n, d, R, x)
    # to account for min/max
    if x[2] == 1.0 || x[2] == -1.0
        0
    else
        u = 3
        (d*R[o, 2] - n*R[u, 2]) / d^2
    end
end

function dh_dθ(o, n, d, R, x)
    #over under numerator denominator R
    u = 3
    θ = x[3]
    l = x[4]
    w = x[5]
    M = [-sin(θ) -cos(θ)
         cos(θ)  -sin(θ)]

    d_high = ([R[o, 1] R[o, 2]] * M * [l,w])
    d_low = ([R[u, 1] R[u, 2]] * M * [l,w])
    d_high = d_high[1]
    d_low = d_low[1]

    (d/2 * d_high - n/2 * d_low) / d^2
end

function dh_dl(o, n, d, R, x)
    u = 3
    θ = x[3]
    d_high = R[o, 1]*cos(θ) + R[o, 2]sin(θ)
    d_low = R[u, 1]*cos(θ) + R[u, 2]*sin(θ)
    (d/2 * d_high - n/2 * d_low) / d^2
end

function dh_dw(o, n, d, R, x)
    u = 3
    θ = x[3]
    d_high = -R[o, 1]*sin(θ) + R[o, 2]cos(θ)
    d_low = -R[u, 1]*sin(θ) + R[u, 2]*cos(θ)
    (d/2 * d_high - n/2 * d_low) / d^2
end

function dh_dh(o, n, d, R, x)
    h = x[6]
    if h == 0
        0
    else
        u = 3
        (d*R[o, 3] - n*R[u, 3]) / d^2
    end
end
