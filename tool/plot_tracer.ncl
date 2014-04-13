load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 100
    time_step = 1
    file_prefix = "tracers.deform.240x120"

    PI = atan(1.0)*4
    RAD = PI/180

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/Debug/"+file_prefix+".%3.3d.nc\n' $i; done")

    f = addfile(fs(0), "r")

    num_time = dimsizes(fs)
    dims = dimsizes(f->c)
    num_tracer = dims(0)
    delete(dims)
    dims = dimsizes(f->s2)
    num_skel2 = dims(1)
    delete(dims)

    wks = gsn_open_wks("pdf", file_prefix+".tracers")

    point_size = 0.002

    res_map = True
    res_map@gsnFrame = False
    res_map@mpCenterLonF = 180
    res_map@mpOutlineOn = False

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

    ; check = ispan(10000, 10240, 10)
    check = (/28565/)

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        ; res_map@gsnPolar = "NH"
        ; res_map@mpMinLatF = 70
        ; res_map@gsnPolar = "SH"
        ; res_map@mpMaxLatF = -70
        ; if (isdefined("check")) then
        ;     res_map@mpProjection = "Satellite"
        ;     res_map@mpCenterLonF = f->c(check(0),0)/RAD
        ;     res_map@mpCenterLatF = f->c(check(0),1)/RAD
        ;     res_map@mpLimitMode = "Angles"
        ;     res_map@mpLeftAngleF = 2.0
        ;     res_map@mpRightAngleF = res_map@mpLeftAngleF
        ;     res_map@mpTopAngleF = res_map@mpLeftAngleF
        ;     res_map@mpBottomAngleF = res_map@mpLeftAngleF
        ; end if
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        do i = 0, num_tracer-1, 1
            c = (/f->c(i,0),f->c(i,1)/)/RAD
            if (isdefined("check") .and. any(check .eq. i)) then
                res_tracer@gsMarkerColor = "blue"
            else
                res_tracer@gsMarkerColor = "black"
            end if
            gsn_polymarker(wks, map, c(0), c(1), res_tracer)
        end do
        ; tracer shape
        if (.not. isdefined("check")) then
            frame(wks)
            continue
        end if
        do k = 0, dimsizes(check)-1
            i = check(k)
            do j = 0, 3
                s = f->s1(i,j,:)/RAD
                gsn_polymarker(wks, map, s(0), s(1), res_s1)
            end do
            do j = 0, num_skel2-1
                s = f->s2(i,j,:)/RAD
                gsn_polymarker(wks, map, s(0), s(1), res_s2)
            end do
        end do
        frame(wks)
    end do

end
