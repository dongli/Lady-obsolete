load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    load "../build/Debug/neighbors.txt"

    PI = atan(1.0)*4
    RAD = PI/180
    p0 = p0/RAD
    ngb = ngb/RAD
    c0 = c0/RAD

    dims = dimsizes(ngb)
    num_ngb = dims(0)

    dims = dimsizes(c0)
    num_c0 = dims(0)

    wks = gsn_open_wks("pdf", "neighbors")

    res_map = True
    res_map@gsnMaximize = True
    res_map@gsnFrame = False
    res_map@mpProjection = "Satellite"
    res_map@mpLimitMode = "Angles"
    res_map@mpCenterLonF = p0(0)
    res_map@mpCenterLatF = p0(1)
    res_map@mpLeftAngleF = 5.0
    res_map@mpRightAngleF = 5.0
    res_map@mpTopAngleF = 5.0
    res_map@mpBottomAngleF = 5.0

    map = gsn_csm_map(wks, res_map)

    res_p0 = True
    res_p0@gsMarkerIndex = 16
    res_p0@gsMarkerSizeF = 0.002

    res_ngb = True
    res_ngb@gsMarkerIndex = 16
    res_ngb@gsMarkerColor = "red"
    res_ngb@gsMarkerSizeF = 0.002

    res_c0 = True
    res_c0@gsMarkerIndex = 16
    res_c0@gsMarkerColor = "blue"
    res_c0@gsMarkerSizeF = 0.002

    gsn_polymarker(wks, map, p0(0), p0(1), res_p0)
    do i = 0, num_ngb-1
        gsn_polymarker(wks, map, ngb(i,0), ngb(i,1), res_ngb)
    end do
    do i = 0, num_c0-1
        gsn_polymarker(wks, map, c0(i,0), c0(i,1), res_c0)
    end do
    ; gsn_polyline(wks, map, c0(:,0), c0(:,1), res_c0)

    frame(wks)

end
