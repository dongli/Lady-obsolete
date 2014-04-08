load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "../build/Debug/neighbors.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    centroid = centroid/RAD
    ngb_cells = ngb_cells/RAD
    ngb_tracers = ngb_tracers/RAD
    shape = shape/RAD

    dims = dimsizes(ngb_cells)
    num_ngb_cell = dims(0)

    dims = dimsizes(ngb_tracers)
    num_ngb_tracer = dims(0)

    dims = dimsizes(shape)
    num_shape = dims(0)

    wks = gsn_open_wks("pdf", "neighbors")

    res_map = True
    res_map@gsnMaximize = False
    res_map@gsnFrame = False
    res_map@mpProjection = "Satellite"
    res_map@mpLimitMode = "Angles"
    res_map@mpCenterLonF = centroid(0)
    res_map@mpCenterLatF = centroid(1)
    res_map@mpLeftAngleF = 5.0
    res_map@mpRightAngleF = 5.0
    res_map@mpTopAngleF = 5.0
    res_map@mpBottomAngleF = 5.0

    map = gsn_csm_map(wks, res_map)

    point_size = 0.008

    res_centroid = True
    res_centroid@gsMarkerIndex = 16
    res_centroid@gsMarkerSizeF = point_size

    res_ngb_cells = True
    res_ngb_cells@gsMarkerIndex = 16
    res_ngb_cells@gsMarkerColor = "red"
    res_ngb_cells@gsMarkerSizeF = point_size

    res_ngb_tracers = True
    res_ngb_tracers@gsMarkerIndex = 5
    res_ngb_tracers@gsMarkerColor = "magenta"
    res_ngb_tracers@gsMarkerSizeF = point_size

    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsMarkerColor = "blue"
    res_shape@gsLineColor = "blue"
    res_shape@gsMarkerSizeF = point_size
    res_shape@gsLineThicknessF = 2

    do i = 0, num_ngb_cell-1
        gsn_polymarker(wks, map, ngb_cells(i,0), ngb_cells(i,1), res_ngb_cells)
    end do
    do i = 0, num_ngb_tracer-1
        gsn_polymarker(wks, map, ngb_tracers(i,0), ngb_tracers(i,1), res_ngb_tracers)
    end do
    ; do i = 0, num_shape-1
    ;     gsn_polymarker(wks, map, shape(i,0), shape(i,1), res_shape)
    ; end do
    gsn_polyline(wks, map, shape(:,0), shape(:,1), res_shape)
    gsn_polymarker(wks, map, centroid(0), centroid(1), res_centroid)

    frame(wks)

end
