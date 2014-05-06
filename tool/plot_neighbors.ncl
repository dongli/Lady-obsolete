load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "../build/Debug/neighbors.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    point_size = 0.005
    satellite_angle = 2.0

    wks = gsn_open_wks("pdf", "neighbors")

    centroid = centroid/RAD

    res_map = True
    res_map@gsnMaximize = False
    res_map@gsnFrame = False
    ; res_map@mpProjection = "Satellite"
    ; res_map@mpLimitMode = "Angles"
    res_map@mpCenterLonF = 180.0
    ; res_map@mpCenterLonF = centroid(0)
    ; res_map@mpCenterLatF = centroid(1)
    ; res_map@mpLeftAngleF = satellite_angle
    ; res_map@mpRightAngleF = satellite_angle
    ; res_map@mpTopAngleF = satellite_angle
    ; res_map@mpBottomAngleF = satellite_angle
    map = gsn_csm_map(wks, res_map)

    ; --------------------------------------------------------------------------
    ; plot neighboring cells
    if (isdefined("ngb_cells")) then
        ngb_cells = ngb_cells/RAD
        dims = dimsizes(ngb_cells)
        num_ngb_cell = dims(0)
        res_ngb_cells = True
        res_ngb_cells@gsMarkerIndex = 16
        res_ngb_cells@gsMarkerColor = "red"
        res_ngb_cells@gsMarkerSizeF = point_size
        do i = 0, num_ngb_cell-1
            gsn_polymarker(wks, map, ngb_cells(i,0), ngb_cells(i,1), res_ngb_cells)
        end do
    end if

    ; --------------------------------------------------------------------------
    ; plot neighboring tracers
    ; if (isdefined("ngb_tracers")) then
    ;     ngb_tracers = ngb_tracers/RAD
    ;     dims = dimsizes(ngb_tracers)
    ;     num_ngb_tracer = dims(0)
    ;     res_ngb_tracers = True
    ;     res_ngb_tracers@gsMarkerIndex = 5
    ;     res_ngb_tracers@gsMarkerColor = "magenta"
    ;     res_ngb_tracers@gsMarkerSizeF = point_size
    ;     do i = 0, num_ngb_tracer-1
    ;         gsn_polymarker(wks, map, ngb_tracers(i,0), ngb_tracers(i,1), res_ngb_tracers)
    ;     end do
    ; end if

    ; --------------------------------------------------------------------------
    ; plot parcel shape
    shape = shape/RAD
    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsLineColor = "black"
    res_shape@gsLineThicknessF = 1
    gsn_polyline(wks, map, shape(:,0), shape(:,1), res_shape)

    ; --------------------------------------------------------------------------
    ; plot centroid
    res_centroid = True
    res_centroid@gsMarkerIndex = 5
    res_centroid@gsMarkerSizeF = point_size
    gsn_polymarker(wks, map, centroid(0), centroid(1), res_centroid)

    frame(wks)

end
