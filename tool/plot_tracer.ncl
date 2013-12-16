load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    RAD = atan(1.0)/45.0

    fs = systemfunc("ls tracers.nc")

    f = addfile(fs(0), "r")

    num_time = dimsizes(fs)
    dims = dimsizes(f->q)
    num_tracer = dims(0)
    delete(dims)

    wks = gsn_open_wks("pdf", "tracers")

    res_map = True
    res_map@gsnFrame = False

    res_tracer = True
    res_tracer@gsMarkerIndex = 16
    res_tracer@gsMarkerSizeF = 0.002

    do j = 0, num_time-1
        system("echo time step: "+j)
        f = addfile(fs(j), "r")
        map = gsn_csm_map(wks, res_map)
        do i = 0, num_tracer-1
            gsn_polymarker(wks, map, f->q(i,0)/RAD, f->q(i,1)/RAD, res_tracer)
        end do
        frame(wks)
    end do

end
