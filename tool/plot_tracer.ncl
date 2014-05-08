load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    file_prefix = "tracers.deform.sc.240x120.split20"
    run_type = "Debug"
    time_stamp = "%4.4d"
    start_time = 130
    end_time = 138
    time_step = 1

    PI = atan(1.0)*4
    RAD = PI/180

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '../build/"+run_type+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

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
    res_map@mpCenterLonF = 180.0
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

    ; check = ispan(7500, 7500, 100)
    ; check = (/7500,12066,28240,28241,46364,46365,70250,70251/)
    check = (/20256,24420,24421/)
    ; colors = (/"red","blue","blue","green","green","magenta","magenta"/)

    s2 = new((/2,num_skel2+1/), double)
    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        if (isdefined("id")) then
            delete(id)
        end if
        if (isdefined("c")) then
            delete(c)
        end if
        id = f->id
        c = f->c/RAD
        ; res_map@gsnPolar = "NH"
        ; res_map@mpMinLatF = 0
        ; res_map@gsnPolar = "SH"
        ; res_map@mpMaxLatF = -45.0
        if (isdefined("check")) then
            i = ind(id .eq. check(0))
            if (ismissing(i)) then
                continue
            end if
            res_map@mpProjection = "Satellite"
            res_map@mpCenterLonF = doubletofloat(c(i,0))
            res_map@mpCenterLatF = doubletofloat(c(i,1))
            res_map@mpLimitMode = "Angles"
            res_map@mpLeftAngleF = 8.0
            res_map@mpRightAngleF = res_map@mpLeftAngleF
            res_map@mpTopAngleF = res_map@mpLeftAngleF
            res_map@mpBottomAngleF = res_map@mpLeftAngleF
        end if
        ; res_map@mpLimitMode = "LatLon"
        ; res_map@mpMinLonF = 240.0
        ; res_map@mpMaxLonF = 330.0
        ; res_map@mpMinLatF = -70.0
        ; res_map@mpMaxLatF = -10.0
        ; res_map@mpGridAndLimbOn = True
        ; res_map@mpGridLineColor = "grey"
        ; res_map@mpGridLonSpacingF = 1.5
        ; res_map@mpGridLatSpacingF = 1.5
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        ; do i = 0, num_tracer-1, 10
        ;     if (isdefined("check") .and. any(check .eq. id(i))) then
        ;         res_tracer@gsMarkerColor = "blue"
        ;     else
        ;         res_tracer@gsMarkerColor = "black"
        ;     end if
        ;     gsn_polymarker(wks, map, c(i,0), c(i,1), res_tracer)
        ; end do
        ; tracer shape
        if (.not. isdefined("check")) then
            frame(wks)
            continue
        end if
        do k = 0, dimsizes(check)-1
            i = ind(id .eq. check(k))
            if (ismissing(i)) then
                continue
            end if
            ; do j = 0, 3
            ;     s1 = f->s1(i,j,:)/RAD
            ;     gsn_polymarker(wks, map, s1(0), s1(1), res_s1)
            ; end do
            ; res_s2@gsMarkerColor = colors(k)
            ; res_s2@gsLineColor = colors(k)
            do j = 0, num_skel2-1
                s2(:,j) = f->s2(i,j,:)/RAD
                ; gsn_polymarker(wks, map, s2(0,j), s2(1,j), res_s2)
            end do
            s2(:,num_skel2) = s2(:,0)
            gsn_polyline(wks, map, s2(0,:), s2(1,:), res_s2)
        end do

        frame(wks)
    end do

end
