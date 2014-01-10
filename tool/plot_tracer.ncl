load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    PI = atan(1.0)*4
    RAD = PI/180

    fs = systemfunc("ls ../Debug/tracers.600.nc")

    f = addfile(fs(0), "r")

    num_time = dimsizes(fs)
    dims = dimsizes(f->q)
    num_tracer = dims(0)
    delete(dims)

    wks = gsn_open_wks("pdf", "tracers")

    res_map = True
    res_map@gsnMaximize = True
    res_map@gsnFrame = False
    ; res_map@gsnPolar = "SH"
    ; res_map@mpMinLatF = 80
    ; res_map@mpMaxLatF = -80

    res_tracer = True
    res_tracer@gsMarkerIndex = 16
    res_tracer@gsMarkerSizeF = 0.002

    res_shape = True
    res_shape@gsMarkerIndex = 16
    res_shape@gsMarkerColor = "red"
    res_shape@gsMarkerSizeF = 0.002
    res_shape@gsLineColor = "red"

    res_skel = True
    res_skel@gsMarkerIndex = 16
    res_skel@gsMarkerSizeF = 0.002
    res_skel@gsMarkerColor = "green"

    ; n = 20
    ; x = new((/n,2/), double)
    ; y = new((/n,2/), double)
    ; dtheta = 2*PI/(n-1)
    ; do i = 0, n-1
    ;     theta = i*dtheta;
    ;     y(i,0) = cos(theta)
    ;     y(i,1) = sin(theta)
    ; end do

    n = 4
    x = new((/n+1,2/), double)
    y = new((/n,2/), double)
    y(0,:) = (/-0.5, 0.0/)
    y(1,:) = (/ 0.0,-0.5/)
    y(2,:) = (/ 0.5, 0.0/)
    y(3,:) = (/ 0.0, 0.5/)

    ; check = ispan(127, 127, 1)

    do l = 0, num_time-1, 5
        system("echo time step: "+l)
        f = addfile(fs(l), "r")
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        do i = 0, num_tracer-1
            q = (/f->q(i,0),f->q(i,1)/)/RAD
            gsn_polymarker(wks, map, q(0), q(1), res_tracer)
        end do
        ; tracer shape
        do i = 0, num_tracer-1
            if (isdefined("check") .and. .not. any(check .eq. i)) then
                continue
            end if
            do k = 0, n-1
               x(k,:) = f->q(i,:)+f->h(i,:,:)#y(k,:)
               gsn_polymarker(wks, map, x(k,0)/RAD, x(k,1)/RAD, res_shape)
            end do
            x(n,:) = x(0,:)
            x = x/RAD
            gsn_polyline(wks, map, x(:,0), x(:,1), res_shape)
        end do
        ; tracer skeleton
        do i = 0, num_tracer-1
            if (isdefined("check") .and. .not. any(check .eq. i)) then
                continue
            end if
            do j = 0, 3
                s = f->skel(i,j,:)/RAD
                gsn_polymarker(wks, map, s(0), s(1), res_skel)
            end do
        end do
        frame(wks)
    end do

end
