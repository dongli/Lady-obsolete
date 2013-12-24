load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    RAD = atan(1.0)/45.0

    fs = systemfunc("ls ../Debug/velocity.*.nc")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", "velocity")

    res = True
    ; res@gsnPolar = "SH"
    ; res@mpMaxLatF = -85.0
    res@vcRefMagnitudeF = 0.1
    res@vcRefLengthF = 0.05
    res@vcGlyphStyle = "CurlyVector"

    do j = 0, num_time-1, 20
        system("echo time step: "+j)
        f = addfile(fs(j), "r")
        plot = gsn_csm_vector_map(wks, f->u(0::2,0::2), f->v(0::2,0::2), res)
    end do

end
