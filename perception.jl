using LinearAlgebra # for the I

struct ObjectState
    x::Float64
    y::Float64
    θ::Float64
    length::Float64
    width::Float64
    height::Float64
end

struct TracksMessage
    timestamp::Float64
    tracks::Dict{Int, ObjectState}
end



function object_tracker(SENSE::Channel, TRACKS::Channel, EMG::Channel, camera_array, road)
    lines = []

    # see PinholeCamera struct line 60 in sensors.jl

    while true
        sleep(0)
        @return_if_told(EMG)
        meas = @fetch_or_continue(SENSE)
        #tracks = TracksMessage(...)
        #TODO your code here
        #@replace(TRACKS, tracks)    
    end
end

function h_state_to_bbox(full_state, camera)
    # full_state = [x y theta l w h v \omega]
    # in sensors.jl:
    # get_camera_meas -> expected_bbox
    #                 -> get_corners
    # get_corners is in movables.jl line 98
# ∈
    # ms stands for moveables
    # you get moveables from simulator_state which is from SIM_ALL in launch files
    # wait you can look at the moveables generated in launch_perception starting from line 24
    # 
    # update_ensor with camera array is the one being used, not the piholecamera one

    (x, y, θ, l, w, h) = full_state[1:6]

    f = camera.focal_len
    R = camera.R
    t = camera.t

    A = [cos(θ)/2  -sin(θ)/2 0
         sin(θ)/2  cos(θ)/2 0
         0         0        1]
    # question: shouldnt the 2nd row sin be negative according ot notes???

    lwh = [l, w, h]
    pos = [x, y, 0]
    points = []
    # + and - lwh to get all 8 points for possible BB coords
    for i in [1, -1]
        for j in [1, -1]
            for k in [0, 1]
                # rf_cam = R*(pos + A*Diagonal([i,j, k])*lwh) + t
                lwh_config = [i, j, k]
                # include as a tuple the corner and also a label for which corner it is
                push!(points, (R*(pos + A*Diagonal(lwh_config)*lwh) + t, lwh_config))
            end
        end
    end

    (left, top, right, bottom) = (Inf, Inf, -Inf, -Inf)
    (left_pt, top_pt, right_pt, bottom_pt) = (0, 0, 0, 0)
    counter = 1
    # look at expected_bbox in sensors.jl, line 117
    for (rf_cam,) in points
        # px = max(min(camera.focal_len * pt[1] / (pt[3] * camera.sx), 1.0), -1.0)
        cam_x = max(min(f*rf_cam[1]/(rf_cam[3]*camera.sx), 1.0), -1.0)
        cam_y = max(min(f*rf_cam[2]/(rf_cam[3]*camera.sy), 1.0), -1.0)
        if cam_x > right
            right = cam_x
            right_pt = counter
        end
        if cam_y < top
            top = cam_y
            top_pt = counter
        end
        if cam_x < left
            left = cam_x
            left_pt = counter
        end
        if cam_y > bottom
            bottom = cam_y
            bottom_pt = counter
        end
        counter += 1
    end

    bbox = [left top right bottom]
    bbox_i = [left_pt top_pt right_pt bottom_pt]

    return (bbox, bbox_i, points) # if just returning h output
end

# function test_J()
function test_J(full_state, camera)
    # full_state = randn(8)
    # camera = 
    (bbox, bbox_i, points) = h_state_to_bbox(full_state, camera)
    jacobian = h_jacobian(bbox_i, points, camera, full_state)

    d = 1e-4

    for i in 1:8
        s = copy(full_state)
        s[i] += d
        (bbox_, _, _) = h_state_to_bbox(s, camera)
        num_jac = (bbox_ - bbox)./d
        
        println("numerical: ")
        println(num_jac)
        println("my calculated:")
        println(jacobian[:,i])
        println()
    end
    println()
    println()
    println()
    println()
    println()
end

function h_jacobian(bbox_i, points, camera, full_state)
    j = zeros(4,8)
    (left, top, right, bottom) = bbox_i
    R = camera.R
    count = 1
    # everything * f/camera.sx or sy
    for i in [left, right]
        (pt, lwh_c) = points[i]
        # x -> 1/3
        n = pt[1]
        d = pt[3]

        config = ones(8)
        config[4:6] = lwh_c
        s = copy(full_state) .* config # shows which corner

        row = zeros(8)
        row[1] = dh_dx(1, n, d, R) # x
        row[2] = dh_dy(1, n, d, R) # y
        row[3] = dh_dθ(1, n, d, R, s)
        row[4] = dh_dl(1, n, d, R, s)
        row[5] = dh_dw(1, n, d, R, s)
        row[6] = dh_dh(1, n, d, R, s)
        row .* camera.focal_len/camera.sx
        j[count, :] = row
        count += 1
    end
    for i in [top, bottom]
        (pt, lwh_c) = points[i]
        # x -> 2/3
        n = pt[2]
        d = pt[3]

        config = ones(8)
        config[4:6] = lwh_c
        s = copy(full_state) .* config

        row = zeros(8)
        row[1] = dh_dx(2, n, d, R) # x
        row[2] = dh_dy(2, n, d, R) # y
        row[3] = dh_dθ(2, n, d, R, s)
        row[4] = dh_dl(2, n, d, R, s)
        row[5] = dh_dw(2, n, d, R, s)
        row[6] = dh_dh(2, n, d, R, s)
        row .* camera.focal_len/camera.sy
        j[count, :] = row
        count += 1
    end

    row2 = copy(j[2, :])
    j[2, :] = j[3, :]
    j[3, :] = row2
    # swaprow!(j,2,3) # so that it is left top right bottom instead of left right top bottom
    return j
end

function dh_dx(o, n, d, R)
    u = 3
    (d*R[o, 1] - n*R[u, 1]) / d^2
end

function dh_dy(o, n, d, R)
    u = 3
    (d*R[o, 2] - n*R[u, 2]) / d^2
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

function obj_state_next(obj_state, cam_meas)
    # x_k_forecast + kalman_gain * (cam_meas - h(x_k_forecast))
    x_f = obj_state_forecast(obj_state, 0.1) # what is time step
end

function obj_state_forecast(obj_state, Δ)
    # dynamics!
    θ = obj_state[3]
    v = obj_state[7]
    ω = obj_state[8]
    obj_state[1:3] +=  [Δ * cos(θ) * v, Δ * cos(θ) * v, Δ * ω]
    return obj_state
end

function J_dynamics_forecast(obj_state, Δ)
    # jacobian of above.. should be 8x8 (dimensions of object state)
end

function kalman_gain()
    # p_forecast * J_h(x_k_forecast).T * (J_h(x_k_forecast) * P_f * J_h(x_k_forecast).T + noise_R)^-1
end

# p is covariance
function p_k()
    # (I - kalman_gain*J_h(x_forecast))*p_forecast
end

function p_forecast(prev_state)
    # J_dynamics(prev_state)*p_prev*J_dynamics(prev_state).T + noise_Q
end
