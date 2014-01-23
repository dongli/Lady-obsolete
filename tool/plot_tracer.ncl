load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    PI = atan(1.0)*4
    RAD = PI/180

    fs = systemfunc("ls ../build/Debug/tracers.300.nc")

    f = addfile(fs(0), "r")

    num_time = dimsizes(fs)
    dims = dimsizes(f->q)
    num_tracer = dims(0)
    delete(dims)

    wks = gsn_open_wks("pdf", "tracers")

    point_size = 0.002

    res_tracer = True
    res_tracer@gsMarkerIndex = 16
    res_tracer@gsMarkerSizeF = point_size

    res_s1 = True
    res_s1@gsMarkerIndex = 16
    res_s1@gsMarkerSizeF = point_size
    res_s1@gsMarkerColor = "green"

    res_s2 = True
    res_s2@gsMarkerIndex = 16
    res_s2@gsMarkerColor = "red"
    res_s2@gsMarkerSizeF = point_size
    res_s2@gsLineColor = "red"

    ; check = ispan(115, 115, 1)
    check = (/10/)

    do l = 0, num_time-1, 2
        system("echo time step: "+l)
        f = addfile(fs(l), "r")
        res_map = True
        res_map@gsnMaximize = True
        res_map@gsnFrame = False
        ; res_map@gsnPolar = "SH"
        ; res_map@mpMaxLatF = -60
        ; if (isdefined("check")) then
        ;     res_map@mpProjection = "Satellite"
        ;     res_map@mpCenterLonF = f->q(check(0),0)/RAD
        ;     res_map@mpCenterLatF = f->q(check(0),1)/RAD
        ;     res_map@mpLimitMode = "Angles"
        ;     res_map@mpLeftAngleF = 8.0
        ;     res_map@mpRightAngleF = res_map@mpLeftAngleF
        ;     res_map@mpTopAngleF = res_map@mpLeftAngleF
        ;     res_map@mpBottomAngleF = res_map@mpLeftAngleF
        ; end if
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        do i = 0, num_tracer-1
            q = (/f->q(i,0),f->q(i,1)/)/RAD
            if (isdefined("check") .and. any(check .eq. i)) then
                res_tracer@gsMarkerColor = "blue"
            else
                res_tracer@gsMarkerColor = "black"
            end if
            gsn_polymarker(wks, map, q(0), q(1), res_tracer)
        end do
        ; tracer shape
        do i = 0, num_tracer-1
            if (isdefined("check") .and. .not. any(check .eq. i)) then
                continue
            end if
            do j = 0, 3
                s = f->s1(i,j,:)/RAD
                gsn_polymarker(wks, map, s(0), s(1), res_s1)
            end do
            x = new((/5,2/), double)
            do j = 0, 3
                x(j,:) = f->s2(i,j,:)/RAD
                gsn_polymarker(wks, map, x(j,0), x(j,1), res_s2)
            end do
            x(4,:) = x(0,:)
            gsn_polyline(wks, map, x(:,0), x(:,1), res_s2)
        end do
        frame(wks)
    end do

end
