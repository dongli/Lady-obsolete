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
    res_map@gsnFrame = False
    ; res_map@gsnPolar = "SH"
    ; res_map@mpMinLatF = 80
    ; res_map@mpMaxLatF = -70

    res_tracer = True
    res_tracer@gsMarkerIndex = 16
    res_tracer@gsMarkerSizeF = 0.002

    res_shape = True
    res_shape@gsLineColor = "red"

    n = 20
    x = new((/n,2/), double)
    y = new((/n,2/), double)
    dtheta = 2*PI/(n-1)
    do i = 0, n-1
        theta = i*dtheta;
        y(i,0) = 0.5*cos(theta)
        y(i,1) = 0.5*sin(theta)
    end do
    do j = 0, num_time-1, 20
        system("echo time step: "+j)
        f = addfile(fs(j), "r")
        map = gsn_csm_map(wks, res_map)
        ; tracer centroid
        do i = 0, num_tracer-1
            q = (/f->q(i,0),f->q(i,1)/)/RAD
            gsn_polymarker(wks, map, q(0), q(1), res_tracer)
        end do
        ; tracer shape
        do i = 0, num_tracer-1
            do k = 0, n-1
               x(k,:) = f->q(i,:)+f->h(i,:,:)#y(k,:)
            end do
            x = x/RAD
            gsn_polyline(wks, map, x(:,0), x(:,1), res_shape)
        end do
        frame(wks)
    end do

end
