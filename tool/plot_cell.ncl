load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    fs = systemfunc("for (( i = 0; i <= 240; ++i )); do printf '../build/Debug/tracers.%3.3d.nc\n' $i; done")

    num_time = dimsizes(fs)

    wks = gsn_open_wks("pdf", "cells")

    gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.0
    res@cnMaxLevelValF = 750.0
    res@cnLevelSpacingF = 50.0
    res@gsnPolar = "SH"
    res@mpMaxLatF = -20

    do l = 0, num_time-1
        system("echo time step: "+l)
        f = addfile(fs(l), "r")
        res@gsnRightString = "Max:"+sprintf("%8.1f", max(f->q0(1:179,:)))
        plot = gsn_csm_contour_map(wks, f->q0(1:179,:), res)
    end do

end
