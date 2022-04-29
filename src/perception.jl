using LinearAlgebra # for the I

struct ObjectState
    x::Float64
    y::Float64
    θ::Float64
    length::Float64
    width::Float64
    height::Float64
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

    while true
        sleep(0)
        @return_if_told(EMG)
        meas = @fetch_or_continue(SENSE)
        
        #tracks = TracksMessage(...)
        #TODO your code here
        #@replace(TRACKS, tracks)    
    end
end

# initialize values for first iteration of Kalman filter
function kalman_init(bb_c1, c1; loop_radius=50.0)
    # c1=camera_array[1]
    f = c1.focal_len
    R = c1.R
    t = c1.t
    sx = c1.sx
    sy = c1.sy
    camera_pos = SVector{3,Float64}(5.0/6*loop_radius, loop_radius, 20.0)
    lookat = SVector{3,Float64}(0, 0, 0)


    #bb_c1 = bb_sense[1]
    obj_estimates = [] # global frame 3d points representing center of each bb in c1

    for bbox in bb_c1
        println(bbox)
        println()
        tl = [bbox.left bbox.top]
        br = [bbox.right bbox.bottom]
        center_pix = midpoint(tl, br)
        center_cf = zeros(Float64, 3,)
        center_cf[3] = euclidean_dist(camera_pos, lookat)
        center_cf[1] = center_pix[1] * (sx/f) * center_cf[3]
        center_cf[2] = center_pix[2] * (sy/f) * center_cf[3]
        println(center_cf)
        println()
        println(t)
        println()

        center_gf = inv(R)*(center_cf - t)
        println(center_gf)
        println()
        x = center_gf[1]
        y = center_gf[2]
        extra_size = rand()
        w = 1.5 + 2.0 * extra_size
        l = 3.0 + 6.0 * extra_size
        h = 2.0 + 1.0 * extra_size
        v = 5.0 + rand()*5.0
        θ = rand() * 2.0 * pi - pi
        angular_vel = v / loop_radius
        push!(obj_estimates, [x y θ l w h v angular_vel])
    end

    println(obj_estimates)
    println()
    println()

    obj_estimates
end

function match_bb(bb_cf, camera_array)
    # transform input bb in global frame
    # match by euclidean distance, output stacked Zk
    # BEFORE: transform everything into global reference frame points
        # x3* = 1.0 for all... thus, x1* = y1/f, x2* = y2/f
        # [x1*, x2*, x3*] --> into [x1, x2, x3] --> create bounding boxes assuming f = 1 
        # y1 = x1/x3, y2 = x2/x3 --> compare these below
    # RETURN the transformed global bounding boxes
    # 1. calculate euclidean distance between TL and BR, add them
    # 2. put into matrix (n x n) where a12 is the distance between the first and second bb
    # note: this is dependent on whether we know which camera is associated with which bb
    # 3. find minimum in each row and col, if multiple match with same box, give priority to the lowest match
    # 4. 

    # get global points of bounding boxes
    c1_bb_global = []
    c2_bb_global = []
    for i in 1:length(bb_cf[1])
        bb1 = camera_to_global(camera_array[1], bb_cf[1][i])
        push!(bb1, c1_bb_global)
    end
    for i in 1:length(bb_cf[2])
        bb2 = camera_to_global(camera_array[2], bb_cf[2][i])
        push!(bb2, c2_bb_global)
    end

    # construct matrix of (euclidean distances)^2
    # stores sum of TL and BR euclidean dist between ith meas from c1 and jth meas from c2
    row = length(c1_bb_global)
    col = length(c2_bb_global)
    dist_matrix = zeros(Float64, row, col)
    min_row = [] # keep (smallest dist, ind) in a row
    min_col = [] # keep (smallest dist, ind) in a col
    for (i, c1_bb) in enumerate(c1_bb_global)
        for (j, c2_bb) in enumerate(c2_bb_global)
            p1_tl = c1_bb[1:2]
            p1_br = c1_bb[3:4]
            p2_tl = c2_bb[1:2]
            p2_br = c2_bb[3:4]
            p1_midpoint = midpoint(p1_tl, p1_br)
            p2_midpoint = midpoint(p2_tl, p2_br)
            sum_dist = euclidean_dist_sq(p1_midpoint, p2_midpoint)
            dist_matrix[i, j] = sum_dist
        end
        # findmin gives (value, index) of min element; index is corresponding bb match
        push!(min_row, findmin(dist_matrix[i]))
    end
    for j in 1:length(c2_bb_global)
        push!(min_col, findmin(dist_matrix[:,j]))
    end

    # add match if both row and col found each other as min
    matches = []
    bb_final = []
    last_ind = max(row, col)
    for i in 1:last_ind
        row_match = min_row[i][2]
        if min_col[row_match][2] == i 
            push!(matches, [i, row_match])
            c1_bb = bb_cf[1][i]
            c2_bb = bb_cf[2][row_match]
            push!(bb_final, [c1_bb c2_bb])
        end
    end

    bb_final
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
