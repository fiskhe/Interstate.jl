using LinearAlgebra # for the I

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

        # for the jacobian, set pt index to invalid number so we know j is zero
        if cam_x == 1.0
            right_pt = 0
        elseif cam_x == -1.0
            left_pt = 0
        end

        if cam_y == 1.0
            bottom_pt = 0
        elseif cam_x == -1.0
            top_pt = 0
        end
        counter += 1
    end

    bbox = [left top right bottom]
    bbox_i = [left_pt top_pt right_pt bottom_pt]

    return (bbox, bbox_i, points) # if just returning h output
end

function test_J(full_state, camera)
    (bbox, bbox_i, points) = h_state_to_bbox(full_state, camera)
    jacobian = h_jacobian_deconstr(bbox_i, points, camera, full_state)

    d = 1e-4
    println("bbox: ")
    println(bbox)

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

function h_jacobian_deconstr(bbox_i, points, camera, full_state)
    f = camera.focal_len
    sx = camera.sx
    sy = camera.sy
    jacobian = zeros(4, 8)
    for (index, i) in enumerate(bbox_i)
        if i != 0
            (pt, lwh_c) = points[i]
            config = ones(8)
            config[4:6] = lwh_c
            s = copy(full_state) .* config # shows which corner
            (th, l, w, h) = s[3:6]

            di_dcorners = zeros(8)
            di_dcorners[i] = 1

            dcorner_dpt = zeros(8, 3)
            if index % 2 == 1 # account for if x or y coordinate
                dcorner_dpt[i, :] = [f/(sx*pt[3]) 0 -f*pt[1]/(sx*pt[3]^2)]
            else
                dcorner_dpt[i, :] = [0 f/(sy*pt[3]) -f*pt[2]/(sy*pt[3]^2)]
            end
            dpt_dtransform = camera.R

            dtrans_dx = [1, 0, 0]
            dtrans_dy = [0, 1, 0]
            dtrans_dth = [-l*sin(th)/2-w*cos(th)/2, l*cos(th)/2-w*sin(th)/2, 0]
            dtrans_dl = [cos(th)/2, sin(th)/2, 0]
            dtrans_dw = [-sin(th)/2, cos(th)/2, 0]
            dtrans_dh = [0, 0, 1]
            if h == 0
                dtrans_dh[3] = 0
            end
            dtr_ds = hcat(dtrans_dx, dtrans_dy, dtrans_dth, dtrans_dl, dtrans_dw, dtrans_dh, zeros(3, 2))

            i_jac = di_dcorners' * dcorner_dpt * dpt_dtransform * dtr_ds
            jacobian[index, :] = i_jac
        end
    end
    jacobian
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

function kalman_gain(prev_state)
    # p_forecast * J_h(x_k_forecast).T * (J_h(x_k_forecast) * P_f * J_h(x_k_forecast).T + noise_R)^-1
    # p_forecast(prev_state) * 
end

# p is covariance
function p_k()
    # (I - kalman_gain*J_h(x_forecast))*p_forecast
end

function p_forecast(prev_state)
    # J_dynamics(prev_state)*p_prev*J_dynamics(prev_state).T + noise_Q
end

