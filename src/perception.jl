using LinearAlgebra # for the I

struct ObjectState
    x::Float64
    y::Float64
    θ::Float64
    length::Float64
    width::Float64
    height::Float64
end

struct ExtendedObjectState
    x::Float64
    y::Float64
    θ::Float64
    length::Float64
    width::Float64
    height::Float64
    velocity::Float64
    ang_velocity::Float64
end

mutable struct KalmanState
    x_k::Vector{Float64}
    P_k::Matrix{Float64}
    meas::Vector{Float64}
    timestamp::Float64
end

function position(o::ObjectState)
    [o.x, o.y]
end

function rear(o::ObjectState)
    -o.length/2.0
end
   
function front(o::ObjectState)
    o.length/2.0
end

function heading(o::ObjectState)
    o.θ
end

struct TracksMessage
    timestamp::Float64
    tracks::Dict{Int, ObjectState}
end


function object_tracker(SENSE::Channel, TRACKS::Channel, EMG::Channel, camera_array, road)
    lines = []

    # kalman filter initialization
    c1 = camera_array[1]
    c2 = camera_array[2]

    first = false
    prev_t = 0
    kalman_states = Vector{KalmanState}()

    while true
        sleep(0)
        @return_if_told(EMG)
        meas = @fetch_or_continue(SENSE)
        if length(meas[1]) == 0 && length(meas[2]) == 0
            continue
        end

        c1_meas = meas[1]
        c2_meas = meas[2]

        if !first
            # initialize vector of all current tracks
            kalman_init(kalman_states, c1_meas, c1)
            prev_t = kalman_states[1].timestamp
            first = true
            continue
        end

        println("before filter: ")
        println(kalman_states)
        println()

        c1_meas = copy(meas[1])
        c2_meas = copy(meas[2])

        for state in kalman_states
            if(length(c1_meas) > 0)
                c1_match = match_bb_to_track(state, c1_meas, c1)
                c1_meas_tmp = c1_meas[c1_match]
                c1_match_meas = [c1_meas_tmp.left, c1_meas_tmp.top, c1_meas_tmp.right, c1_meas_tmp.bottom]
                cur_t = c1_meas_tmp.time
                deleteat!(c1_meas, c1_match)
            else
                c1_match_meas = [0, 0, 0, 0]
            end
            if(length(c2_meas) > 0)
                c2_match = match_bb_to_track(state, c2_meas, c2)
                c2_meas_tmp = c2_meas[c2_match]
                c2_match_meas = [c2_meas_tmp.left, c2_meas_tmp.top, c2_meas_tmp.right, c2_meas_tmp.bottom]
                cur_t = c2_meas_tmp.time
                deleteat!(c2_meas, c2_match)
            else 
                c2_match_meas = [0, 0, 0, 0]
            end
            state.meas = vcat(c1_match_meas, c2_match_meas)
            delta_t = cur_t - prev_t
            state.timestamp = cur_t
            kalman_filter(state, camera_array, delta_t)
            prev_t = cur_t
        end

        # add new tracks for leftover bounding boxes in camera1
        if length(c1_meas) != 0
            kalman_init(kalman_states, c1_meas, c1)
        end

        println("after filter: ")
        println(kalman_states)
        println()
        
        #tracks = TracksMessage(...)
        #TODO your code here
        #@replace(TRACKS, tracks)    
    end
end

function kalman_filter(prev_kalman_state, cams, delta_t)
   
    x_prev = prev_kalman_state.x_k
    P_prev = prev_kalman_state.P_k

    x_k_forecast = obj_state_forecast(x_prev, 0.1) # 0.1 time step for now

    c1 = cams[1]
    c2 = cams[2]
    (bb1, bbox_i1, points1) = h_state_to_bbox(x_k_forecast, c1)
    (bb2, bbox_i2, points2) = h_state_to_bbox(x_k_forecast, c2)
    bb_forecast = vcat(bb1', bb2')

    J_h1= h_jacobian_deconstr(bbox_i1, points1, c1, x_prev)
    J_h2= h_jacobian_deconstr(bbox_i2, points2, c2, x_prev)
    J_h = vcat(J_h1, J_h2)

    P_k_forecast = p_forecast(x_prev, P_prev, delta_t)
    kalman = kalman_gain(P_k_forecast, x_k_forecast, J_h)
    P_k = p_k(kalman, x_k_forecast, P_k_forecast, J_h)
    x_k = obj_state_next(x_k_forecast, kalman, prev_kalman_state.meas, bb_forecast)

    prev_kalman_state.x_k = vec(x_k)
    prev_kalman_state.P_k = P_k

end

# initialize values for first iteration of Kalman filter with camera1
function kalman_init(kalman_states, bb_c1, c1; loop_radius=50.0)
    # c1=camera_array[1]
    f = c1.focal_len
    R = c1.R
    t = c1.t
    sx = c1.sx
    sy = c1.sy
    camera_pos = SVector{3,Float64}(5.0/6*loop_radius, loop_radius, 20.0)
    lookat = SVector{3,Float64}(0, 0, 0)

    for (i, bbox) in enumerate(bb_c1)
        tl = [bbox.left bbox.top]
        br = [bbox.right bbox.bottom]
        center_pix = midpoint(tl, br)
        center_cf = zeros(Float64, 3,)
        center_cf[3] = euclidean_dist(camera_pos, lookat)
        center_cf[1] = center_pix[1] * (sx/f) * center_cf[3]
        center_cf[2] = center_pix[2] * (sy/f) * center_cf[3]

        center_gf = inv(R)*(center_cf - t)
        x = center_gf[1]
        y = center_gf[2]
        extra_size = rand()
        w = 1.5 + 2.0 * extra_size
        l = 3.0 + 6.0 * extra_size
        h = 2.0 + 1.0 * extra_size
        v = 5.0 + rand()*5.0
        θ = rand() * 2.0 * pi - pi
        angular_vel = v / loop_radius
        x_0 = [x, y, θ, l, w, h, v, angular_vel]
        variances = [25, 25, 25, 25, 4, 1, 25, 0.04] # all overestimations
        P_0 = diagm(variances)

        curr_state = KalmanState(x_0, P_0, [bbox.left, bbox.top, bbox.right, bbox.bottom], bbox.time)
        push!(kalman_states, curr_state)
    end
end

function match_bb_to_track(kalman_state, bb_list, camera)
    (kalman_bbox, _, _) = h_state_to_bbox(kalman_state.x_k, camera)
    euclid_dists = []
    for (i, bb) in enumerate(bb_list)
        tl_kalman = [kalman_bbox[1] kalman_bbox[2]]
        br_kalman = [kalman_bbox[3] kalman_bbox[4]]
        kalman_bbox_mid = midpoint(tl_kalman, br_kalman)
        tl_bbox = [bb.left bb.top]
        br_bbox = [bb.right bb.bottom]
        bb_mid = midpoint(tl_bbox, br_bbox)
        push!(euclid_dists, sqrt((kalman_bbox_mid[1]-bb_mid[1])^2 + (kalman_bbox_mid[2]-bb_mid[2])^2))
    end
    best_match = findmin(euclid_dists)
    return best_match[2] #index of best match bb
end


# converts bounding boxes from camera frame to global frame 2d pixels (assumes f=1)
function camera_to_global(camera, bb_cf)
    f = camera.focal_len
    R = camera.R
    t = camera.t

    # get 3d coordinates in camera reference frame (TL and BR)
    x1_tl = bb_cf[1]/f 
    x2_tl = bb_cf[2]/f
    x3_tl = 1.0
    cf_3d_tl = SVector{3,Float64}(x1_tl, x2_tl, x3_tl)
    x1_br = bb_cf[3]/f
    x2_br = bb_cf[4]/f 
    x3_br = 1.0
    cf_3d_br = SVector{3,Float64}(x1_br, x2_br, x3_br)

    # get 3d coordinates in global reference frame
    gf_3d_tl = inv(R)*(cf_3d_tl - t)
    gf_3d_br = inv(R)*(cf_3d_br - t)

    # convert to bb pixels with f = 1
    y1_tl = gf_3d_tl[1]/gf_3d_tl[3]
    y2_tl = gf_3d_tl[2]/gf_3d_tl[3]
    y1_br = gf_3d_br[1]/gf_3d_br[3]
    y2_br = gf_3d_br[2]/gf_3d_br[3]
    return [y1_tl y2_tl y1_br y2_br]
end

# returns the square of euclidean distance between 2 points
function euclidean_dist(p1, p2)
    return sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2 + (p1[3]-p2[3])^2)
end

# returns the center between 2 points
function midpoint(p1, p2)
    return [(p1[1]+p2[1])/2 (p1[2]+p2[2])/2]
end
